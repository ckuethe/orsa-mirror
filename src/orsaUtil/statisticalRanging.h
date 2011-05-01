#ifndef _ORSA_UTIL_STATISTICAL_RANGING_H_
#define _ORSA_UTIL_STATISTICAL_RANGING_H_

#include <orsaUtil/adaptiveInterval.h>
#include <orsaUtil/adaptiveMonteCarlo.h>
#include <orsaUtil/multimin.h>

#include <orsa/statistic.h>

#include <orsaSolarSystem/observation.h>
#include <orsaSolarSystem/orbit.h>

namespace orsaUtil {

    class SR_AuxiliaryData : public osg::Referenced {
    protected:
        virtual ~SR_AuxiliaryData() { }
    public:
        orsaSolarSystem::OpticalObservationVector allOpticalObs;
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
                            const double & residualProbability, // 1-confidenceLevel
                            const double & initialThresholdLevel,
                            const double & targetThresholdLevel,
                            const size_t & targetSamples,
                            const SR_AuxiliaryData * auxiliaryData) :
            orsaUtil::AdaptiveInterval<SR_ElementData> (min,max,residualProbability,initialThresholdLevel,targetThresholdLevel,targetSamples),
            aux(auxiliaryData)
            { }
    protected:
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
                    const orsa::Time t = aux->allOpticalObs[j]->epoch;
                    orsaSolarSystem::OrbitWithEpoch tmpOrbit = e.data->O_s2an_g;
                    tmpOrbit.M = fmod(e.data->O_s2an_g.M + orsa::twopi()*(t-e.data->O_s2an_g.epoch).get_d()/e.data->O_s2an_g.period(),orsa::twopi());
                    tmpOrbit.relativePosVel(rOrbit,vOrbit);
                    e.data->R_o2an[j] = (rOrbit+aux->R_s[j]-aux->R_o[j]);
                    e.data->V_o2an[j] = (vOrbit+aux->V_s[j]-aux->V_o[j]);
                }
#warning need to do it here...
                // update position
                if (!e.position.isLocked()) {
#warning should set if already set? even if not locked?
                    e.position = (*e.data->R_o2an[j]).length();
                    e.position.lock();
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
    
    typedef std::vector< osg::ref_ptr<SR_AdaptiveInterval> > SR_AdaptiveIntervalVector;
    
    class SR_AdaptiveMonteCarlo : public orsaUtil::AdaptiveMonteCarlo<SR_AdaptiveIntervalVector> {
    public:
        SR_AdaptiveMonteCarlo(const SR_AuxiliaryData * auxiliaryData) :
            orsaUtil::AdaptiveMonteCarlo<SR_AdaptiveIntervalVector>(),
            aux(auxiliaryData),
            mov(new orsaUtil::MultiminOrbitalVelocity) 
            { }
    protected:
        ~SR_AdaptiveMonteCarlo() { }
    protected:
        osg::ref_ptr<const SR_AuxiliaryData> aux;
        osg::ref_ptr<orsaUtil::MultiminOrbitalVelocity> mov;
        
        // two different integers between 0 and N-1, with z1 < z2
        /* static bool z1z2(size_t & z1, size_t & z2, const size_t & N) {
           if (N<2) return false;
           // don't use ___z1 and ___z2 in code below, but z1 and z2 
           size_t ___z1 = orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform_int(N);
           size_t ___z2;
           do { ___z2 = orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform_int(N); }
           while (___z1 == ___z2);
           // save them as z1 < z2 to ensure that vectors...[z1] are set before vectors...[z2]
           z1 = std::min(___z1,___z2);
           z2 = std::max(___z1,___z2);
           return true;
           }
        */
        // for this one, we don't force z1 < z2, so it can be z1 > z2
    public:
        static bool z1z2(size_t & z1, size_t & z2, const size_t & N) {
            if (N<2) return false;
            z1 = orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform_int(N);
            do { z2 = orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform_int(N); }
            while (z1 == z2);
            return true;
        }
    public:
        bool sample(ElementTypeVector & ev) const {
            std::vector<orsa::Vector> R_s2an;
            R_s2an.resize(aux->vecSize);
            size_t z1=0, z2=0; // set to 0 just to make compiler happy...
            bool bound;
            double dx,dy;
#warning a local maxIter should be added, as the do/while loop could be an infinite loop in certain situations
            do {
                bound=true;
                // note: z1 and z2 are not ordered, i.e. can be z2 < z1
                z1z2(z1,z2,aux->vecSize);
                // ORSA_DEBUG("z1: %i  z2: %i",z1,z2);
                // for (unsigned int j=0; j<aux->vecSize; ++j) {
                // if ((j!=z1) && (j!=z2)) continue;

                unsigned int j;
                j = z1;
                // if (j==z1)
                {
                    orsa::GlobalRNG::instance()->rng()->gsl_ran_dir_2d(&dx,&dy);
                    const orsa::Vector u_o2an_j = (aux->u_o2a[j]+orsa::GlobalRNG::instance()->rng()->gsl_ran_gaussian(aux->vecRMS[j]*orsa::arcsecToRad())*(dx*aux->xS_o2a[j]+dy*aux->yS_o2a[j])).normalized();
                    const orsa::Vector R_o2an_j = intervalVector[j]->sample()*u_o2an_j;
                    const orsa::Vector R_an_j = aux->R_o[j] + R_o2an_j;
                    R_s2an[j] = R_an_j - aux->R_s[j];
                    
                }
                // if (j==z2)
                j = z2;
                {
                    // enforce bound orbit on z2
                    const double escapeVelocity = sqrt(2*orsaSolarSystem::Data::GMSun()/R_s2an[z1].length());
                    const double dt = fabs((aux->allOpticalObs[z2]->epoch-aux->allOpticalObs[z1]->epoch).get_d());
                    const double maxDistance = escapeVelocity*dt;
                    orsa::GlobalRNG::instance()->rng()->gsl_ran_dir_2d(&dx,&dy);
                    const orsa::Vector u_o2an_j = (aux->u_o2a[j]+orsa::GlobalRNG::instance()->rng()->gsl_ran_gaussian(aux->vecRMS[j]*orsa::arcsecToRad())*(dx*aux->xS_o2a[j]+dy*aux->yS_o2a[j])).normalized();
                    const orsa::Vector diffVector = R_s2an[z1]-aux->R_o[z2]+aux->R_s[z2];
                    const double rangeCenter = std::max(0.0,u_o2an_j*diffVector); // range of min distance to R...z1
                    const double minDistance = (u_o2an_j*rangeCenter-diffVector).length();
                    if (minDistance > maxDistance) {
                        // no possible bound solutions here
                        bound=false;
                        // continue;
                    } else {
                        const double rangeDelta  = sqrt(orsa::square(maxDistance)-orsa::square(minDistance));
                        const orsa::Vector R_o2an_j = std::max(0.0,rangeCenter+rangeDelta*(2*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform()-1))*u_o2an_j;
                        const orsa::Vector R_an_j = aux->R_o[j] + R_o2an_j;
                        R_s2an[j] = R_an_j - aux->R_s[j];
                    }
                }
                // if (j==z2) break; // done // can't do this because z1 and z2 are not ordered...
                // }
            } while (bound==false);

            // ORSA_DEBUG("bound: %i", bound);
            
            const orsa::Vector V_s2an_z1 = mov->getOrbitalVelocity(R_s2an[z1],aux->allOpticalObs[z1]->epoch,
                                                                   R_s2an[z2],aux->allOpticalObs[z2]->epoch,
#warning replace this with GM with something like aux->GM_s;
                                                                   orsaSolarSystem::Data::GMSun());
            
            {
                orsaSolarSystem::OrbitWithEpoch tmpOrbit;
                tmpOrbit.compute(R_s2an[z1],
                                 V_s2an_z1,
                                 orsaSolarSystem::Data::GMSun());
                tmpOrbit.epoch = aux->allOpticalObs[z1]->epoch;
                const orsaSolarSystem::OrbitWithEpoch O_s2an_g = tmpOrbit;
                
                for (unsigned int k=0; k<aux->vecSize; ++k) {
                    ev[k].data->resize(aux->vecSize);
                    ev[k].data->O_s2an_g = O_s2an_g;
                }
            }
            
            for (unsigned int k=0; k<aux->vecSize; ++k) {
                intervalVector[k]->updateLevel(ev[k]);
            }
            
            return true;
        } 
    };
    
    bool statisticalRanging(SR_AdaptiveIntervalVector & vec,
                            const double & initialThresholdLevel,
                            const double & targetThresholdLevel,
                            const double & minAdaptiveRange,
                            const double & maxAdaptiveRange,
                            const double & intervalResidualProbability,
                            const size_t & targetSamples,
                            const size_t & maxIter,
                            const SR_AuxiliaryData * aux);
    
}; // namespace orsaUtil

#endif // _ORSA_UTIL_STATISTICAL_RANGING_H_
