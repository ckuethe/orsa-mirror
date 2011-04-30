#ifndef _ORSA_UTIL_STATISTICAL_RANGING_H_
#define _ORSA_UTIL_STATISTICAL_RANGING_H_

#include <orsaUtil/adaptiveInterval.h>
#include <orsaUtil/adaptiveMonteCarlo.h>

#include <orsa/statistic.h>

#include <orsaSolarSystem/observation.h>
#include <orsaSolarSystem/orbit.h>

namespace orsaUtil {

    class SR_AuxiliaryData : public osg::Referenced {
    protected:
        virtual ~SR_AuxiliaryData() { }
    public:
        std::vector< orsa::Cache<orsa::Vector> > R_s, V_s, R_o, V_o;
        std::vector< orsa::Cache<orsa::Vector> > u_o2a; // unit vectors, from Dawn to Satellite
        std::vector< orsa::Cache<orsa::Vector> > xS_o2a, yS_o2a; // unit vectors, orthogonal to u_d2s, to model astrometric accuacy
        std::vector<double> vecRMSnominal; // same index as observations
        std::vector<double> vecRMS; // same index as observations
        orsa::Cache<size_t> vecSize;
    };

    class SR_ElementData : public osg::Referenced {
    public:
        void resize(const size_t & vecSize) {
            vec_residual.resize(vecSize);
            R_o2an.resize(vecSize);
            V_o2an.resize(vecSize);
        }
    protected:
        ~SR_ElementData() { }
    public:
        orsaSolarSystem::OrbitWithEpoch O_s2an_g;
    public:
        mutable std::vector< orsa::Cache<double> > vec_residual;
        mutable orsa::Cache<double> RMS;
    public:
        mutable std::vector< orsa::Cache<orsa::Vector> > R_o2an;
        mutable std::vector< orsa::Cache<orsa::Vector> > V_o2an;
    };
    
    typedef orsaUtil::AdaptiveIntervalElement<SR_ElementData> SR_Element;
    
    class SR_AdaptiveInterval : public orsaUtil::AdaptiveInterval<SR_ElementData> {
    public:
        SR_AdaptiveInterval(const double & min,
                            const double & max,
                            const double & confidenceLevel,
                            const double & initialThresholdLevel,
                            const double & targetThresholdLevel,
                            const size_t & targetSamples,
                            const orsaSolarSystem::OpticalObservationVector & allOpticalObs_,
                            const SR_AuxiliaryData * auxiliaryData) :
            orsaUtil::AdaptiveInterval<SR_ElementData> (min,max,confidenceLevel,initialThresholdLevel,targetThresholdLevel,targetSamples),
            allOpticalObs(allOpticalObs_),
            aux(auxiliaryData)
            { }
    protected:
        const orsaSolarSystem::OpticalObservationVector & allOpticalObs;
        osg::ref_ptr<const SR_AuxiliaryData> aux;
    public:
        void updateLevel(const AdaptiveIntervalElementType & e) const {
            if (e.level.isSet()) return;
            osg::ref_ptr< orsa::Statistic<double> > stat_residual =
                new orsa::Statistic<double>;
            double chisq=0.0; // chisq
            orsa::Vector rOrbit, vOrbit;
            for (unsigned int j=0; j<aux->vecSize; ++j) {
                if ( (!e.data->R_o2an[j].isSet()) ||
                     (!e.data->V_o2an[j].isSet()) ) {
                    const orsa::Time t = allOpticalObs[j]->epoch;
                    orsaSolarSystem::OrbitWithEpoch tmpOrbit = e.data->O_s2an_g;
                    tmpOrbit.M = fmod(e.data->O_s2an_g.M + orsa::twopi()*(t-e.data->O_s2an_g.epoch).get_d()/e.data->O_s2an_g.period(),orsa::twopi());
                    tmpOrbit.relativePosVel(rOrbit,vOrbit);
                    e.data->R_o2an[j] = (rOrbit+aux->R_s[j]-aux->R_o[j]);
                    e.data->V_o2an[j] = (vOrbit+aux->V_s[j]-aux->V_o[j]);
                }
                if (!e.data->vec_residual[j].isSet()) {
                    const orsa::Vector u_o2an_j = (*e.data->R_o2an[j]).normalized();
                    const double residual = acos(std::min(1.0,u_o2an_j*aux->u_o2a[j]))*orsa::radToArcsec();
                    e.data->vec_residual[j] = residual;
                }
                stat_residual->insert(e.data->vec_residual[j]);
#warning use vecRMS or vecRMSnominal? if using nominal (that never changes) then updating level is not necessary...
                chisq+=orsa::square(e.data->vec_residual[j]/aux->vecRMS[j]);
            }
            e.level = chisq;
            e.data->RMS = stat_residual->RMS();
        }    
    };
    
    bool statisticalRanging(const orsaSolarSystem::OpticalObservationVector & allOpticalObs,
                            const SR_AuxiliaryData * aux);
    
    
}; // namespace orsaUtil

#endif // _ORSA_UTIL_STATISTICAL_RANGING_H_
