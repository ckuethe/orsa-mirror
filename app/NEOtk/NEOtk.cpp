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
#include <orsaUtil/multimin.h>

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
                         const orsaSolarSystem::OpticalObservationVector & allOpticalObs_) :
        orsaUtil::AdaptiveInterval<AdaptiveIntervalTemplateType> (min,max,confidenceLevel,initialThresholdLevel,targetThresholdLevel),
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
            chisq+=orsa::square(e.data->vec_residual[j]/vecRMS[j]);
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
        
        const double minAdaptiveRange = orsa::FromUnits(100.0,orsa::Unit::KM); // important: not zero...
        const double maxAdaptiveRange = orsa::FromUnits(100.0,orsa::Unit::AU);
        
        const double intervalResidualProbability = 1.0e-10; // "1-confidence level" for this interval, different from the chisq-level
        
        const size_t targetSamples = 100;
        
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
        
        std::vector< orsa::Cache<double> > vecMinRange;
        vecMinRange.resize(allOpticalObs.size());
        std::vector< orsa::Cache<double> > vecMaxRange;
        vecMaxRange.resize(allOpticalObs.size());
        {
            // test: determine the maximum range for each observation
            osg::ref_ptr<orsaUtil::MultiminMinMaxRange> mmr = new orsaUtil::MultiminMinMaxRange;
            
            // run on ALL j and k because range search is asymmetric (determines at t=tj)
            for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
                    if (j==k) continue;
                    orsa::Cache<double> minRange, maxRange;
                    mmr->getMinMaxRange(minRange, maxRange,
                                        minAdaptiveRange,
                                        maxAdaptiveRange,
                                        R_s[j],R_o[j],u_o2a[j],vecRMS[j],allOpticalObs[j]->epoch,
                                        R_s[k],R_o[k],u_o2a[k],vecRMS[k],allOpticalObs[k]->epoch,
                                        orsaSolarSystem::Data::GMSun(),
                                        3.0);
#warning nominal max error here... 3.0? 2.5?
                    
                    // maxRange is only for first set of arguments == [j]
                    // vecMaxRange[j].setIfLarger(maxRange);
                    // if (maxRange > 0.0) {
                    if (minRange.isSet()) {
                        if (minRange > 0.0) {
                            vecMinRange[j].setIfSmaller(minRange);
                        }
                    }
                    if (maxRange.isSet()) {
                        if (maxRange > 0.0) {
                            vecMaxRange[j].setIfLarger(maxRange);
                        }
                    }

                    /* if (maxRange.isSet()) {
                       ORSA_DEBUG("max range distance for obs [%02i] and [%02i] = %8.3g [AU]",
                       j,
                       k,
                       orsa::FromUnits(maxRange,orsa::Unit::AU,-1));
                       }
                    */
                    
                }
                // ORSA_DEBUG("set: %i",vecMaxRange[j].isSet());
                if (vecMinRange[j].isSet()) {
                    ORSA_DEBUG("min range distance for obs [%02i] = %8.3f [AU]",
                               j,
                               orsa::FromUnits(vecMinRange[j],orsa::Unit::AU,-1));
                } else {
#warning is it correct to force it here?
                    vecMinRange[j] = minAdaptiveRange; 
                }
                if (vecMaxRange[j].isSet()) {
                    ORSA_DEBUG("max range distance for obs [%02i] = %8.3f [AU]",
                               j,
                               orsa::FromUnits(vecMaxRange[j],orsa::Unit::AU,-1));
                } else {
#warning is it correct to force it here?
                    vecMaxRange[j] = maxAdaptiveRange; 
                }
            }           
        }
        
        std::vector< osg::ref_ptr<AdaptiveIntervalType> > range_99;
        range_99.resize(allOpticalObs.size());
        for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
            range_99[k] = new AdaptiveIntervalType(std::max(minAdaptiveRange,(*vecMinRange[k])),
                                                   std::min((*vecMaxRange[k]),maxAdaptiveRange),
                                                   intervalResidualProbability, // "1-confidence level" for this interval, different from the chisq-level
                                                   chisq_99,
                                                   chisq_99,
                                                   allOpticalObs);
        }
        
        std::vector< osg::ref_ptr<AdaptiveIntervalType> > range_XY;
        range_XY.resize(allOpticalObs.size());
        for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
            range_XY[k] = new AdaptiveIntervalType(std::max(minAdaptiveRange,(*vecMinRange[k])),
                                                   std::min((*vecMaxRange[k]),maxAdaptiveRange),
                                                   intervalResidualProbability, // "1-confidence level" for this interval, different from the chisq-level
#warning this level should depend on the success, i.e. start with low level and increase it if no points get inserted...
                                                   1000*chisq_99,
                                                   chisq_99, // this should be the lowest possible level we work with
                                                   allOpticalObs);
        }
        
        std::vector< osg::ref_ptr<AdaptiveIntervalType> > activeRange = range_XY;
        
        std::vector<size_t> activeRangeOldSize;
        activeRangeOldSize.resize(vecSize);
        for (unsigned int k=0; k<vecSize; ++k) {
            activeRangeOldSize[k] = 0;
        }
        
        std::vector<size_t> activeRangeOldSizeDBG;
        activeRangeOldSizeDBG.resize(vecSize);
        for (unsigned int k=0; k<vecSize; ++k) {
            activeRangeOldSizeDBG[k] = 0;
        }
        
        osg::ref_ptr<orsaUtil::MultiminOrbitalVelocity> mov = new orsaUtil::MultiminOrbitalVelocity;
        unsigned int iter=0;
        // test counts
#warning need to improve this! (one for each value of chisq confidence level?)
        unsigned int ct_tot=0;
        unsigned int ct_NEO=0;
        unsigned int old_ct_tot=0;
        
        orsa::Cache<double> minRMS;
        while (iter<100000000) {
            ++iter;            
            
            std::vector<orsa::Vector> R_s2an; R_s2an.resize(allOpticalObs.size());
            // std::vector<orsa::Vector> V_s2an; V_s2an.resize(allOpticalObs.size());
            
            // don't use ___z1 and ___z2 in code below, but z1 and z2 
            unsigned int ___z1 = orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform_int(allOpticalObs.size());
            unsigned int ___z2;
            do { ___z2 = orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform_int(allOpticalObs.size()); }
            while (___z1 == ___z2);
            // save them as z1 < z2 to ensure that vectors...[z1] are set before vectors...[z2]
            const unsigned int z1 = std::min(___z1,___z2);
            const unsigned int z2 = std::max(___z1,___z2);
            {
                bool bound=true;
                double dx,dy;
                for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                    if ((j!=z1) && (j!=z2)) continue;
                    if (j==z1) {
                        orsa::GlobalRNG::instance()->rng()->gsl_ran_dir_2d(&dx,&dy);
                        const orsa::Vector u_o2an_j = (u_o2a[j]+orsa::GlobalRNG::instance()->rng()->gsl_ran_gaussian(vecRMS[j]*orsa::arcsecToRad())*(dx*xS_o2a[j]+dy*yS_o2a[j])).normalized();
                        const orsa::Vector R_o2an_j = activeRange[j]->sample()*u_o2an_j;
                        const orsa::Vector R_an_j = R_o[j] + R_o2an_j;
                        R_s2an[j] = R_an_j - R_s[j];
                        
                    }
                    if (j==z2) {
                        // enforce bound orbit on z2
                        const double escapeVelocity = sqrt(2*orsaSolarSystem::Data::GMSun()/R_s2an[z1].length());
                        const double dt = (allOpticalObs[z2]->epoch-allOpticalObs[z1]->epoch).get_d();
                        const double maxDistance = escapeVelocity*dt;
                        orsa::GlobalRNG::instance()->rng()->gsl_ran_dir_2d(&dx,&dy);
                        const orsa::Vector u_o2an_j = (u_o2a[j]+orsa::GlobalRNG::instance()->rng()->gsl_ran_gaussian(vecRMS[j]*orsa::arcsecToRad())*(dx*xS_o2a[j]+dy*yS_o2a[j])).normalized();
                        const orsa::Vector diffVector = R_s2an[z1]-R_o[z2]+R_s[z2];
                        const double rangeCenter = std::max(0.0,u_o2an_j*diffVector); // range of min distance to R...z1
                        const double minDistance = (u_o2an_j*rangeCenter-diffVector).length();
                        if (minDistance > maxDistance) {
                            // no possible bound solutions here
                            bound=false;
                            continue;
                        }
                        const double rangeDelta  = sqrt(orsa::square(maxDistance)-orsa::square(minDistance));
                        const orsa::Vector R_o2an_j = std::max(0.0,rangeCenter+rangeDelta*(2*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform()-1))*u_o2an_j;
                        const orsa::Vector R_an_j = R_o[j] + R_o2an_j;
                        R_s2an[j] = R_an_j - R_s[j];
                    }
                    if (j==z2) break; // done
                }
                if (!bound) continue;
            }
            
            const orsa::Vector V_s2an_z1 = mov->getOrbitalVelocity(R_s2an[z1],allOpticalObs[z1]->epoch,
                                                                   R_s2an[z2],allOpticalObs[z2]->epoch,
                                                                   orsaSolarSystem::Data::GMSun());
            
            {
                orsaSolarSystem::OrbitWithEpoch tmpOrbit;
                tmpOrbit.compute(R_s2an[z1],
                                 V_s2an_z1,
                                 orsaSolarSystem::Data::GMSun());
                tmpOrbit.epoch = allOpticalObs[z1]->epoch;
                const orsaSolarSystem::OrbitWithEpoch O_s2an_g = tmpOrbit;
                
                AdaptiveIntervalElementType e;
                e.data->resize(vecSize);
                e.data->O_s2an_g = O_s2an_g;
                for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
                    activeRange[k]->updateLevel(e);
                    e.position = (*e.data->R_o2an[k]).length();
                    activeRange[k]->insert(e);
                }
            }
            
#warning should NOT be [0] but the one with the smallest number of entries... because for that one, all entries are in common with all other activeRange
            AdaptiveIntervalType::DataType::DataType::const_iterator it = activeRange[0]->getData()->getData().begin();
            while (it != activeRange[0]->getData()->getData().end()) {
                
                // perform this check only if enought points are already available, to speed up things...
#warning should replace this threshold (16) with a check that enough points are available with lower RMS value
                if (iter%100==0) {
                    if (activeRange[0]->size() >= 5) {
                        AdaptiveIntervalElementType newMinRMSEntry;
                        bool newMinRMS=false;
                        
                        if (!minRMS.isSet()) {
                            minRMS = (*it).data->RMS;
                            newMinRMSEntry=(*it);
                            newMinRMS=true;
                        } else {
                            if ((*it).data->RMS < minRMS) {
                                minRMS = (*it).data->RMS;
                                newMinRMSEntry=(*it);
                                newMinRMS=true;
                            }
                        }
                        
                        if (newMinRMS) {
                            ORSA_DEBUG("new minRMS = %g",(*minRMS));
#warning remember to reset all other counters around the code...
                            for (unsigned int k=0; k<allOpticalObs.size(); ++k) {

                                // not here, just a little bit later...
                                // activeRange[k]->refresh();
                                
                                // reset counters...
                                ct_tot=0;
                                ct_NEO=0;
                                old_ct_tot=0;
                                
                                //
                                
#warning RE-INCLUDE THIS!
                                {
                                    const double old_vecRMS_k = vecRMS[k];
                                    if (newMinRMSEntry.data->vec_residual[k] >= vecRMSnominal[k]) {
#warning MAXIMUM RMS value set here, with min(...,3.0)...
                                        // vecRMS[k] = std::min((*newMinRMSEntry.data->vec_residual[k]),2.5);
                                        vecRMS[k] = std::min((*newMinRMSEntry.data->vec_residual[k]),2*vecRMSnominal[k]);
                                    } else {
                                        vecRMS[k] = vecRMSnominal[k];
                                    }
                                    
                                    if (vecRMS[k] < old_vecRMS_k) {
                                        ORSA_DEBUG("RMS[%02i] = %5.2f (lowered)",k,vecRMS[k]);
                                    } else if (vecRMS[k] > old_vecRMS_k) {
                                        ORSA_DEBUG("RMS[%02i] = %5.2f (raised)",k,vecRMS[k]);
                                    } else {
                                        ORSA_DEBUG("RMS[%02i] = %5.2f",k,vecRMS[k]);
                                    }
                                }
                            }
                            
                            for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
                                activeRange[k]->refresh();
                            }
                            
                            // continue;
                            break;
                        }
                    }
                }
                
                if (iter%100==0) {

                    // debug
                    
                    bool printSize=false;
                    for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
                        if (activeRange[k]->size()!=activeRangeOldSizeDBG[k]) {
                            printSize=true;
                            break;
                        }
                    }
                    if (printSize) {
                        char line_printSize[1024*1024];
                        char tmpStr[1024];
                        line_printSize[0] = '\0';
                        for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
                            sprintf(tmpStr,"%i %.1f ",activeRange[k]->size(),activeRange[k]->getThreshold());
                            strcat(line_printSize,tmpStr);
                        }
                        ORSA_DEBUG("activeRanges: %s",line_printSize);
                    }
                    // update
                    for (unsigned int k=0; k<vecSize; ++k) {
                        activeRangeOldSizeDBG[k]=activeRange[k]->size();
                    }
                }
                
#warning which level to use?
                // if (stat_residual->RMS() < successThreshold) {
                // if ((*it).level < chisq_99) {
                if ((*it).level < activeRange[0]->getThreshold()) {
                    // stats
                    ++ct_tot;
                    if ((*it).data->O_s2an_g.a*(1.0-(*it).data->O_s2an_g.e) < FromUnits(1.3,orsa::Unit::AU)) {
                        ++ct_NEO;
                    }
                }
                
                ++it;
            }
            
#warning should test this only if there has been an insert...       
            // another test: can decreast threshold level?
            if (iter%100==0) {
                for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
#warning should test if the adaptive interval actually shrinks when reducing level and at the same time reducing the number of points...
                    if (activeRange[k]->size()==activeRangeOldSize[k]) continue;
                    activeRangeOldSize[k]=activeRange[k]->size();
                    if (activeRange[k]->getThreshold() == activeRange[k]->getTargetThreshold()) continue;
                    const double currentSampleRange = activeRange[k]->maxSample()-activeRange[k]->minSample();
#warning parameters...
                    // const size_t minSize = 50;
                    const double thresholdIncreaseFactor = 1.2;
                    // if (activeRange[k]->size()<minSize) continue;
                    double testThreshold = activeRange[k]->getTargetThreshold();
                    double testSampleRange = currentSampleRange;
                    orsa::Cache<double> testMin, testMax;
                    size_t testSize = 0;
                    while (1) {
                        if (testThreshold>=activeRange[k]->getThreshold()) break;
                        testSize = 0;
                        {
                            AdaptiveIntervalType::DataType::DataType::const_iterator it =
                                activeRange[k]->getData()->getData().begin();
                            while (it != activeRange[k]->getData()->getData().end()) {
                                if ((*it).level < testThreshold) {
                                    testMin.setIfSmaller((*it).position);
                                    testMax.setIfLarger((*it).position);
                                    ++testSize;
                                }
                                ++it;
                            }
                        }
                        // test if new interval would be smaller than current one
                        if (testMin.isSet() && testMax.isSet()) {
                            // should use an actual AdaptiveInterval for this?
                            
                            // same as 2*delta in
                            // must use same residual probability as actual intervals
                            if (testSize >= 2) {
                                testSampleRange = (testMax-testMin)/(pow(intervalResidualProbability,1.0/testSize));
                                if (testSampleRange < currentSampleRange) {
                                    break; // done
                                }
                            }
                        }
                        // exception break: have enough points at target threshold
                        if ( (testThreshold==activeRange[k]->getTargetThreshold()) && 
                             (testSize>=targetSamples) ) {
                            ORSA_DEBUG("QUICK EXIT for activeRange[%i] ****************************",k);
                            break;
                        }
                        testThreshold *= thresholdIncreaseFactor;
                    }
                    bool reduceThreshold=false;
                    if ( (testSampleRange < currentSampleRange) && 
                         (testSize>=2) &&
                         (testThreshold<activeRange[k]->getThreshold()) ) {
                        reduceThreshold=true;
                    }
                    if ( (testThreshold==activeRange[k]->getTargetThreshold()) &&
                         (testSize>=targetSamples) ) {
                        reduceThreshold=true;
                    }
                    if (reduceThreshold) {
                        ORSA_DEBUG("reducing threshold for range %3i: %8.2f -> %8.2f   size: %3i -> %3i   range: %7.3f -> %7.3f [AU]",
                                   k,
                                   activeRange[k]->getThreshold(),
                                   testThreshold,
                                   activeRange[k]->size(),
                                   testSize,
                                   orsa::FromUnits(currentSampleRange,orsa::Unit::AU,-1),
                                   orsa::FromUnits(testSampleRange,orsa::Unit::AU,-1));
                        activeRange[k]->updateThresholdLevel(testThreshold);
                    }
                }
            }
            
            {
                bool canStop=true;
                for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
                    if ( (activeRange[k]->size() >= targetSamples) &&
                         (activeRange[k]->getThreshold()==activeRange[k]->getTargetThreshold()) ) {
                        // keep going
                    } else {
                        canStop=false;
                    }
                }
                if (canStop) break;
            }
        }
        
        ORSA_DEBUG("iter: %i",iter);
        
#warning should NOT be [0] but the one with the smallest number of entries... because for that one, all entries are in common with all other activeRange
        
        AdaptiveIntervalType::DataType::DataType::const_iterator it = activeRange[0]->getData()->getData().begin();
        while (it != activeRange[0]->getData()->getData().end()) {
            
            // main output
#warning which level to use?
            // if ((*it).level < chisq_99) {
            if ((*it).level < activeRange[0]->getThreshold()) {
                
                const orsaSolarSystem::OrbitWithEpoch O_s2an_g = (*it).data->O_s2an_g;
                
                // a string to write all the single residuals
                char line_eachResidual[1024*1024];
                {
                    char tmpStr[1024];
                    line_eachResidual[0] = '\0';
                    for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                        if ((*it).data->vec_residual[j].isSet()) {
                            sprintf(tmpStr,"%6.2f ",(*(*it).data->vec_residual[j]));
                            strcat(line_eachResidual,tmpStr);
                        }
                    }
                }
                
                // this is relative to the observer
                char line_eachDistance[1024*1024];
                {
                    char tmpStr[1024];
                    line_eachDistance[0] = '\0';
                    for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                        orsa::Cache<orsa::Vector> & cv = (*it).data->R_o2an[j];
                        if (cv.isSet()) {
                            sprintf(tmpStr,"%6.3f ",orsa::FromUnits((*cv).length(),orsa::Unit::AU,-1));
                            strcat(line_eachDistance,tmpStr);
                        }
                    }
                }
                
                char line_eachVelocity[1024*1024];
                {
                    char tmpStr[1024];
                    line_eachVelocity[0] = '\0';
                    for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                        orsa::Cache<orsa::Vector> & cv = (*it).data->V_o2an[j];
                        if (cv.isSet()) {
                            sprintf(tmpStr,"%4.1f ",orsa::FromUnits(orsa::FromUnits((*cv).length(),orsa::Unit::KM,-1),orsa::Unit::SECOND));
                            strcat(line_eachVelocity,tmpStr);
                        }
                    }
                }
                
                // a string to write all the single RMS values used
                char line_eachRMS[1024*1024];
                {
                    char tmpStr[1024];
                    line_eachRMS[0] = '\0';
                    for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                        sprintf(tmpStr,"%6.3f ",vecRMS[j]);
                        strcat(line_eachRMS,tmpStr);
                    }
                }
                
                // absolute magnitude estimate
                double H;
                {
                    osg::ref_ptr< orsa::Statistic<double> > stat_H = new orsa::Statistic<double>;
                    for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                        orsa::Cache<orsa::Vector> & cv = (*it).data->R_o2an[j];
                        if (cv.isSet()) {
                            if (allOpticalObs[j]->mag.isSet()) {
#warning see the correction here to the reported magnitude, depending on the filter 
                                /* stat_H->insert(absoluteMagnitude(allOpticalObs[j]->mag+orsaSolarSystem::MPC_band_correction(allOpticalObs[j]->band),
                                   0.15,
                                   acos((R_o[j]-R_an[j]).normalized()*(R_s[j]-R_an[j]).normalized()),
                                   R_o2an[j].length(),
                                   (R_an[j]-R_s[j]).length()));
                                */
                                const orsa::Vector R_an = cv + R_o[j];
                                stat_H->insert(absoluteMagnitude(allOpticalObs[j]->mag+orsaSolarSystem::MPC_band_correction(allOpticalObs[j]->band),
                                                                 0.15,
                                                                 acos((R_o[j]-R_an).normalized()*(R_s[j]-R_an).normalized()),
                                                                 (*cv).length(),
                                                                 (R_an-R_s[j]).length()));
                            }
                        }
                    }
                    // ORSA_DEBUG("entries: %Zi",stat_H->entries().get_mpz_t());
                    H = stat_H->average();
                }
                
                // Earth MOID
                double EarthMOID;
                {
                    orsa::Orbit EarthOrbit;
                    EarthOrbit.compute(earth.get(),sun.get(),bg.get(),O_s2an_g.epoch);
                    double M1, M2;
                    orsa::MOID(EarthMOID,
                               M1,
                               M2,
                               EarthOrbit,
                               O_s2an_g,
                               16,
                               1.0e-6);
                }
                
                
                ORSA_DEBUG("SAMPLE[minRMS=%.3f]: %6.3f %7.3f %9.3f %8.6f %7.3f %5.2f %6.4f %s %s %s %s",
                           (*minRMS),
                           (*(*it).data->RMS),
                           (*(*it).level),
                           orsa::FromUnits(O_s2an_g.a,orsa::Unit::AU,-1),
                           O_s2an_g.e,
                           O_s2an_g.i*orsa::radToDeg(),
                           H,
                           orsa::FromUnits(EarthMOID,orsa::Unit::AU,-1),
                           line_eachResidual,
                           line_eachDistance,
                           line_eachVelocity,
                           line_eachRMS);
                
            } else {
                // ORSA_DEBUG("above threshold: %g",stat_residual->RMS());
            }
            
            
            if (1) {
                // debug output
                // if (iter%10==0) {
                if (ct_tot!=old_ct_tot) {
                    if (ct_tot!=0) {
                        ORSA_DEBUG("prob. NEO: %7.3f (%i/%i)",(double)ct_NEO/(double)ct_tot,ct_NEO,ct_tot);
                    }
                    old_ct_tot=ct_tot;
                }
            }
            
            ++it;
        }
    }
    
    ORSA_DEBUG("done.");
    
    return 0;
}
