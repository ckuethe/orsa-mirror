#include <orsaSPICE/spice.h>
#include <orsaSPICE/spiceBodyTranslationalCallback.h>

#include <orsaSolarSystem/attitude.h>
#include <orsaSolarSystem/data.h>
#include <orsaSolarSystem/datetime.h>
#include <orsaSolarSystem/obleq.h>
#include <orsaSolarSystem/orbit.h>
#include <orsaSolarSystem/print.h>

#include <orsaUtil/observatory.h>

#include <orsaInputOutput/MPC_observations.h>
#include <orsaInputOutput/MPC_obscode.h>
#include <orsaInputOutput/RWO.h>

#include <orsa/bodygroup.h>
#include <orsa/multimin.h>
#include <orsa/orbit.h>
#include <orsa/print.h>
#include <orsa/integrator_radau.h>
#include <orsa/util.h>
#include <orsa/statistic.h>

#include <gsl/gsl_cdf.h>

#include <algorithm>

#include <QApplication>

#include <orsaUtil/adaptiveInterval.h>
#include <orsaUtil/adaptiveMonteCarlo.h>
#include <orsaUtil/multimin.h>
#include <orsaUtil/statisticalRanging.h>

using namespace orsa;

#warning improve this with something better than a bunch of global vars...
// some global vars (bad!)
// all resized and initialized early in main, before main loop
// s = sun, o = observer, a = asteroid
std::vector< orsa::Cache<orsa::Vector> > R_s, V_s, R_o, V_o;
std::vector< orsa::Cache<orsa::Vector> > u_o2a; // unit vectors, from Dawn to Satellite
std::vector< orsa::Cache<orsa::Vector> > xS_o2a, yS_o2a; // unit vectors, orthogonal to u_d2s, to model astrometric accuacy

std::vector<double> vecRMSnominal; // same index as observations
std::vector<double> vecRMS; // same index as observations
//
orsa::Cache<size_t> vecSize;

class Entry : public osg::Referenced {
#warning use global var for vecSize?...
public:
    void resize(const size_t & vecSize) {
        vec_residual.resize(vecSize);
        R_o2an.resize(vecSize);
        V_o2an.resize(vecSize);
    }
protected:
    ~Entry() { }
public:
    orsaSolarSystem::OrbitWithEpoch O_s2an_g;
public:
    mutable std::vector< orsa::Cache<double> > vec_residual;
    mutable orsa::Cache<double> RMS;
public:
    mutable std::vector< orsa::Cache<orsa::Vector> > R_o2an;
    mutable std::vector< orsa::Cache<orsa::Vector> > V_o2an;
};

typedef Entry AdaptiveIntervalTemplateType;

typedef orsaUtil::AdaptiveIntervalElement<AdaptiveIntervalTemplateType> AdaptiveIntervalElementType;

class AdaptiveIntervalType : public orsaUtil::AdaptiveInterval<AdaptiveIntervalTemplateType> {
public:
    AdaptiveIntervalType(const double & min,
                         const double & max,
                         const double & confidenceLevel,
                         const double & initialThresholdLevel,
                         const double & targetThresholdLevel,
                         const size_t & targetSamples,
                         const orsaSolarSystem::OpticalObservationVector & allOpticalObs_) :
        orsaUtil::AdaptiveInterval<AdaptiveIntervalTemplateType> (min,max,confidenceLevel,initialThresholdLevel,targetThresholdLevel,targetSamples),
        allOpticalObs(allOpticalObs_)
        { }
protected:
    const orsaSolarSystem::OpticalObservationVector & allOpticalObs;
public:
    void updateLevel(const AdaptiveIntervalElementType & e) const {
        if (e.level.isSet()) return;
        osg::ref_ptr< orsa::Statistic<double> > stat_residual =
            new orsa::Statistic<double>;
        double chisq=0.0; // chisq
        orsa::Vector rOrbit, vOrbit;
        for (unsigned int j=0; j<vecSize; ++j) {
            if ( (!e.data->R_o2an[j].isSet()) ||
                 (!e.data->V_o2an[j].isSet()) ) {
                const orsa::Time t = allOpticalObs[j]->epoch;
                orsaSolarSystem::OrbitWithEpoch tmpOrbit = e.data->O_s2an_g;
                tmpOrbit.M = fmod(e.data->O_s2an_g.M + twopi()*(t-e.data->O_s2an_g.epoch).get_d()/e.data->O_s2an_g.period(),orsa::twopi());
                tmpOrbit.relativePosVel(rOrbit,vOrbit);
                e.data->R_o2an[j] = (rOrbit+R_s[j]-R_o[j]);
                e.data->V_o2an[j] = (vOrbit+V_s[j]-V_o[j]);
            }
            if (!e.data->vec_residual[j].isSet()) {  const orsa::Vector u_o2an_j = (*e.data->R_o2an[j]).normalized();
                const double residual = acos(std::min(1.0,u_o2an_j*u_o2a[j]))*orsa::radToArcsec();
                e.data->vec_residual[j] = residual;
            }
            stat_residual->insert(e.data->vec_residual[j]);

#warning divide by vecRMS or vecRMSnominal?
            // chisq+=orsa::square(e.data->vec_residual[j]/vecRMS[j]);
            chisq+=orsa::square(e.data->vec_residual[j]/vecRMSnominal[j]);
            
        }
        e.level = chisq;
        e.data->RMS = stat_residual->RMS();
    }    
};


int main(int argc, char **argv) {
    
    // orsaQt::Debug::instance()->initTimer();
    //
    orsa::Debug::instance()->initTimer();
    
    if (argc != 2) {
        ORSA_DEBUG("Usage: %s <observations-file>",argv[0]);
        exit(0);
    }
    
    orsa::Debug::instance()->initTimer();
  
    orsaSPICE::SPICE::instance()->setDefaultObserver("SSB");
    
    // orsaSPICE::SPICE::instance()->loadKernel("de405.bsp");
    orsaSPICE::SPICE::instance()->loadKernel("de421.bsp");
    
    // read the file obsRMS.dat...
    std::map< std::string, double > obsRMS; // nominal RMS of residuals, in arcsec
    {
        FILE * fp = fopen("obsRMS.dat","r");
        if (!fp) {
            ORSA_DEBUG("problem: cannot open file [obsRMS.dat]");
            exit(0);
        } else {
            char line[1024];
            char obsCode[1024];
            double RMS;
            while (fgets(line,1024,fp)) {
                sscanf(line,"%s %lf",obsCode,&RMS);
                obsRMS[obsCode] = RMS;
            }        
            fclose(fp);
        }
    }
        
    osg::ref_ptr<BodyGroup> bg = new BodyGroup;
    
    osg::ref_ptr<Body> sun = new Body;
    {
        sun->setName("SUN");
        orsaSPICE::SpiceBodyTranslationalCallback * sbtc = 
            new orsaSPICE::SpiceBodyTranslationalCallback(sun->getName());
        orsa::IBPS ibps;
        ibps.inertial = new PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MSun());
        ibps.translational = sbtc;
        sun->setInitialConditions(ibps);
        bg->addBody(sun.get());
    }
    
    osg::ref_ptr<orsa::Body> earth = new orsa::Body;
    {
        earth->setName("EARTH");
        earth->isLightSource = false;
        orsaSPICE::SpiceBodyTranslationalCallback * sbtc = 
            new orsaSPICE::SpiceBodyTranslationalCallback(earth->getName());
        orsa::IBPS ibps;
        ibps.inertial = new orsa::PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MEarth());
        ibps.translational = sbtc;
        earth->setInitialConditions(ibps);
        bg->addBody(earth.get());
    }
    
    osg::ref_ptr<orsa::Body> moon = new orsa::Body;
    {
        moon->setName("MOON");
        moon->isLightSource = false;
        orsaSPICE::SpiceBodyTranslationalCallback * sbtc = 
            new orsaSPICE::SpiceBodyTranslationalCallback(moon->getName());
        orsa::IBPS ibps;
        ibps.inertial = new orsa::PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MMoon());
        ibps.translational = sbtc;
        moon->setInitialConditions(ibps);
        bg->addBody(moon.get());
    }
    
    osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obsCodeFile =
        new orsaInputOutput::MPCObsCodeFile;
    obsCodeFile->setFileName("ObsCodes.html");
    obsCodeFile->read();
    
    osg::ref_ptr<orsaUtil::StandardObservatoryPositionCallback> obsPosCB =
        new orsaUtil::StandardObservatoryPositionCallback(obsCodeFile.get());

    orsaSolarSystem::OpticalObservationVector allOpticalObs;
    //
    {
        osg::ref_ptr< orsaInputOutput::InputFile < orsaInputOutput::CompressedFile,
            std::vector< osg::ref_ptr<orsaSolarSystem::Observation> > > > obsFile;
        //
        {
            osg::ref_ptr<orsaInputOutput::MPCObservationsFile> MPCObsFile = 
                new orsaInputOutput::MPCObservationsFile;
            MPCObsFile->setFileName(argv[1]);
            MPCObsFile->read();
            
            osg::ref_ptr<orsaInputOutput::RWOFile> RWOFile = 
                new orsaInputOutput::RWOFile;
            RWOFile->setFileName(argv[1]);
            RWOFile->read();
            
            const unsigned int maxObs = std::max(MPCObsFile->_data.size(),
                                                 RWOFile->_data.size());
            
            if (maxObs == 0) {
                ORSA_DEBUG("could not read observations file");
                exit (0);
            } else {
                if (maxObs == MPCObsFile->_data.size()) {
                    obsFile = MPCObsFile.get();
                } else if (maxObs == RWOFile->_data.size()) {
                    obsFile = RWOFile.get();
                }
            }
        }
        
        // copy file data to local container
        // observationVector = obsFile->_data;
        
        // copy only optical obs to local container
        for (unsigned int k=0; k<obsFile->_data.size(); ++k) {
            orsaSolarSystem::OpticalObservation * opticalObservation =  
                dynamic_cast<orsaSolarSystem::OpticalObservation *>(obsFile->_data[k].get());
            if (!opticalObservation) {
                ORSA_DEBUG("observation is not an OpticalObservation");
            }
            allOpticalObs.push_back(opticalObservation);
        }
        
        ORSA_DEBUG("total optical observations read: %i",allOpticalObs.size());
        
        // global
        vecSize = allOpticalObs.size();
        
        // assign RMS values
        vecRMSnominal.resize(allOpticalObs.size());
        vecRMS.resize(allOpticalObs.size());
        for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
            double RMS = obsRMS[allOpticalObs[k]->obsCode];
            if (RMS==0.0) {
                // if not set (that is, equal to 0.0), use default
                ORSA_DEBUG("cannot find nominal accuracy for observatory code [%s], please update file obsRMS.dat",(*allOpticalObs[k]->obsCode).c_str());
#warning default RMS=? input parameter?
                RMS=0.9;
            }
            //
            /* allOpticalObs[k]->sigma_ra  = RMS*arcsecToRad();
               allOpticalObs[k]->sigma_dec = RMS*arcsecToRad();
            */
            //
            vecRMSnominal[k] = RMS; // in arcsec
#warning higher initial value of RMS?
            // vecRMS[k]        = 2.50; // in arcsec
            vecRMS[k]        = RMS; // in arcsec
        }
    }
    
    const double chisq_50  = gsl_cdf_chisq_Pinv(0.50,allOpticalObs.size());
    const double chisq_90  = gsl_cdf_chisq_Pinv(0.90,allOpticalObs.size());
    const double chisq_95  = gsl_cdf_chisq_Pinv(0.95,allOpticalObs.size());
    const double chisq_99  = gsl_cdf_chisq_Pinv(0.99,allOpticalObs.size());
    //
    ORSA_DEBUG("chisq 50\%: %.2f  90\%: %.2f  95\%: %.2f  99\%: %.2f",
               chisq_50,chisq_90,chisq_95,chisq_99);
    
    if (1) {
        
        /***** INPUT *****/
        
        // const double minAdaptiveRange = orsa::FromUnits(100.0,orsa::Unit::KM); // important: not zero...
        const double minAdaptiveRange = orsa::FromUnits(  0.0,orsa::Unit::AU); // important: not zero...
        const double maxAdaptiveRange = orsa::FromUnits(100.0,orsa::Unit::AU);
        
        const double intervalResidualProbability = 1.0e-10; // "1-confidence level" for this interval, different from the chisq-level
        
        const size_t targetSamples = 1000;
        
        // const unsigned int minAdaptiveSize = 2; // 2

        // using GlobalRNG
        // const int randomSeed = getpid(); // 717901;
        
        /***** END INPUT *****/
        
        // ORSA_DEBUG("randomSeed: %i",randomSeed);
        // osg::ref_ptr<orsa::RNG> rng = new orsa::RNG(randomSeed);
        
        if (allOpticalObs.size() < 2) {
            ORSA_DEBUG("problem: not enough observations");
            exit(0);
        }
        
        // observation times must be sorted
        for (unsigned int j=1; j<allOpticalObs.size(); ++j) {
            if (allOpticalObs[j]->epoch < allOpticalObs[j-1]->epoch) {
                ORSA_DEBUG("problem: observation times not sorted");
                exit(0);
            }
        }
        
        // s = sun, o = observer, a = asteroid
        // std::vector<orsa::Vector> R_s, V_s, R_o, V_o;
        R_s.resize(allOpticalObs.size());
        V_s.resize(allOpticalObs.size());
        R_o.resize(allOpticalObs.size());
        V_o.resize(allOpticalObs.size());
        
        // std::vector<orsa::Vector> u_o2a; // unit vectors, from Dawn to Satellite
        u_o2a.resize(allOpticalObs.size());
        
        
        // std::vector<orsa::Vector> xS_o2a, yS_o2a; // unit vectors, orthogonal to u_d2s, to model astrometric accuacy
        xS_o2a.resize(allOpticalObs.size());
        yS_o2a.resize(allOpticalObs.size());
        
        for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
            
            const orsa::Time t = allOpticalObs[j]->epoch;

            orsa::Vector tmpR, tmpV;
            
            if (!bg->getInterpolatedPosVel(tmpR,tmpV,sun.get(),t)) {
                ORSA_DEBUG("problems... t:");
                orsa::print(t);
            }
            R_s[j] = tmpR;
            V_s[j] = tmpV;
            
            if (!obsPosCB->getPosVel(tmpR,tmpV,allOpticalObs[j].get())) { 
                ORSA_DEBUG("problems... t:");
                orsa::print(t);
            }
            R_o[j] = tmpR;
            V_o[j] = tmpV;
                
            u_o2a[j] = observationDirection(allOpticalObs[j]);

#warning rotate along RA, Dec?
            xS_o2a[j] = orsa::externalProduct(orsa::Vector(0,0,1),u_o2a[j]).normalized();
            yS_o2a[j] = orsa::externalProduct(u_o2a[j],xS_o2a[j]).normalized();
        } 


        {
// test new code
            osg::ref_ptr<orsaUtil::SR_AuxiliaryData> aux = new orsaUtil::SR_AuxiliaryData;
            aux->allOpticalObs = allOpticalObs;
            aux->R_s = R_s;
            aux->R_o = R_o;
            aux->V_s = V_s;
            aux->V_o = V_o;
            aux->u_o2a = u_o2a;
            aux->xS_o2a = xS_o2a;
            aux->yS_o2a = yS_o2a;
            aux->vecRMSnominal = vecRMSnominal;
            aux->vecRMS = vecRMS;
            aux->vecSize = allOpticalObs.size();
            orsaUtil::SR_AdaptiveIntervalVector vec;
            orsaUtil::statisticalRanging(vec,
                                         100*chisq_99,
                                         chisq_99,
                                         minAdaptiveRange,
                                         maxAdaptiveRange,
                                         intervalResidualProbability,
                                         targetSamples,
                                         1000000, // maxIter
                                         aux.get());

            {
                size_t countALL=0, countNEO=0;
#warning only vec[0] ?      
                orsaUtil::SR_AdaptiveInterval::DataType::DataType::const_iterator it = vec[0]->getData()->getData().begin();
                while (it != vec[0]->getData()->getData().end()) {
                    const orsaSolarSystem::OrbitWithEpoch & o = (*it).data->O_s2an_g;
                    const double q = o.a*(1.0-o.e);
                    if (q < orsa::FromUnits(1.3,orsa::Unit::AU)) {
                        ++countNEO;
                    }
                    ++countALL;
                    ++it;
                }
                ORSA_DEBUG("NEO RATING: %.1f\% (%i/%i)",100*(double)countNEO/(double)countALL,countNEO,countALL);
            }
        }
    }
    
    ORSA_DEBUG("done.");
    
    return 0;
}

