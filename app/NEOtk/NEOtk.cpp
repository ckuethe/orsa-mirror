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

#include <orsaOSG/viz.h>
#include <orsaOSG/FindNamedNodeVisitor.h>
#include <orsaOSG/Track.h>
#include <orsaOSG/DepthPartitionNode.h>
#include <orsaOSG/DistanceAccumulator.h>

#include <osgGA/NodeTrackerManipulator>

#include <gsl/gsl_cdf.h>

#include "vestaViz.h"

#include <algorithm>

#include <QApplication>

#include <orsaQt/debug.h>

#include <orsaUtil/adaptiveInterval.h>
#include <orsaUtil/multimin.h>

using namespace orsa;

// main purpose: disable the rotation after the release of the mouse button
class VestaNodeTrackerManipulator : public osgGA::NodeTrackerManipulator {
public:
    bool handle(const osgGA::GUIEventAdapter & ea,
                osgGA::GUIActionAdapter      & us) {
        if (ea.getEventType() == osgGA::GUIEventAdapter::RELEASE) {
            osg::ref_ptr<osgGA::GUIEventAdapter> mod_ea = new osgGA::GUIEventAdapter(ea);
            mod_ea->setButtonMask(99);
            return osgGA::NodeTrackerManipulator::handle((*mod_ea.get()),us);
        } else {
            return osgGA::NodeTrackerManipulator::handle(ea,us);
        }
    }
};

//

#warning improve this with something better than a bunch of global vars...
// some global vars (bad!)
// all resized and initialized early in main, before main loop
// s = sun, o = observer, a = asteroid
std::vector<orsa::Vector> R_s, V_s, R_o, V_o;
std::vector<orsa::Vector> u_o2a; // unit vectors, from Dawn to Satellite
std::vector<orsa::Vector> xS_o2a, yS_o2a; // unit vectors, orthogonal to u_d2s, to model astrometric accuacy
std::vector<orsa::Vector> u_o2an; // unit vectors, from Dawn to Satellite, with "noise" (astrometric uncertainty)
std::vector<orsa::Vector> R_o2an, V_o2an; // real length vectors, from Dawn to Satellite, with "noise" (astrometric uncertainty)
std::vector<orsa::Vector> R_an, V_an; //
std::vector<orsa::Vector> R_s2an;
std::vector<orsa::Vector> V_s2an;
//
std::vector<double> vecRMSnominal; // same index as observations
std::vector<double> vecRMS; // same index as observations


class Entry {
public:
    orsaSolarSystem::OrbitWithEpoch O_s2an_g;
    // std::vector<double> vec_residual;
    mutable orsa::Cache<double> RMS;
};

typedef Entry AdaptiveIntervalTemplateType;

typedef orsaUtil::AdaptiveIntervalElement<AdaptiveIntervalTemplateType> AdaptiveIntervalElementType;

class AdaptiveIntervalType : public orsaUtil::AdaptiveInterval<AdaptiveIntervalTemplateType> {
public:
    AdaptiveIntervalType(const double & min,
                         const double & max,
                         const double & confidenceLevel,
                         const double & thresholdLevel,
                         const    int & randomSeed,
                         const orsaSolarSystem::OpticalObservationVector & allOpticalObs_) :
        orsaUtil::AdaptiveInterval<AdaptiveIntervalTemplateType> (min,max,confidenceLevel,thresholdLevel,randomSeed),
        allOpticalObs(allOpticalObs_)
        { }
protected:
    const orsaSolarSystem::OpticalObservationVector & allOpticalObs;
protected:
    void updateLevel(const orsaUtil::AdaptiveIntervalElement<Entry> & e) {
        
        if (e.level.isSet()) return;
        
        // compute astrometric offset
        std::vector<double> vec_residual;
        vec_residual.resize(allOpticalObs.size());
        //
        osg::ref_ptr< orsa::Statistic<double> > stat_residual =
            new orsa::Statistic<double>;
        stat_residual->reset();
        double chisq=0.0; // chisq
        orsa::Vector rOrbit, vOrbit;
        for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
            const orsa::Time t = allOpticalObs[j]->epoch;
            orsaSolarSystem::OrbitWithEpoch tmpOrbit = e.O_s2an_g;
            tmpOrbit.M = fmod(e.O_s2an_g.M + twopi()*(t-e.O_s2an_g.epoch).get_d()/e.O_s2an_g.period(),orsa::twopi());
            tmpOrbit.relativePosVel(rOrbit,vOrbit);
            rOrbit += R_s[j];
            vOrbit += V_s[j];
            // debug: print R_sn[j] before changing it
            /* ORSA_DEBUG("distance: %g [km]",
               orsa::FromUnits((R_sn[j]-rOrbit).length(),orsa::Unit::KM,-1));
            */
            // note how some vectors are getting redefined
            R_an[j] = rOrbit;
            V_an[j] = vOrbit;
            R_o2an[j] = (R_an[j]-R_o[j]);
            V_o2an[j] = (V_an[j]-V_o[j]);
            u_o2an[j] = (R_an[j]-R_o[j]).normalized();
            const double residual = acos(std::min(1.0,u_o2an[j]*u_o2a[j]))*orsa::radToArcsec();
            vec_residual[j] = residual;
            stat_residual->insert(residual);
            chisq+=orsa::square(residual/vecRMS[j]);
            // ORSA_DEBUG("OFFSET[%i]: %5.1f",j,residual);
        }
        e.level = chisq;
        e.RMS = stat_residual->RMS();
        // ORSA_DEBUG("chisq: %g  RMS: %g",(*e.level),(*e.RMS));
    }    
};


int main(int argc, char **argv) {
    
    // QApplication app(argc, argv);
    
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
    // test
    /* ORSA_DEBUG("RMS[703] = %g",obsRMS["703"]);
       ORSA_DEBUG("RMS[XYZ] = %g",obsRMS["XYZ"]);
    */
    //
    // std::vector<double> vecRMSnominal; // same index as observations
    // std::vector<double> vecRMS; // same index as observations
    
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
        
        // assign RMS values
        vecRMSnominal.resize(allOpticalObs.size());
        vecRMS.resize(allOpticalObs.size());
        for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
            double RMS = obsRMS[allOpticalObs[k]->obsCode];
            if (RMS==0.0) {
                // if not set (that is, equal to 0.0), use default
                ORSA_DEBUG("cannot find nominal accuracy for observatory code [%s], please update file obsRMS.dat",(*allOpticalObs[k]->obsCode).c_str());
                // #warning default RMS=?   (use automatic, floating RMS?)
                RMS=0.9;
            }
            //
            /* allOpticalObs[k]->sigma_ra  = RMS*arcsecToRad();
               allOpticalObs[k]->sigma_dec = RMS*arcsecToRad();
            */
            //
            vecRMSnominal[k] = RMS; // in arcsec
#warning higher initial value of RMS?
            vecRMS[k]        = 2.50; // in arcsec
        }
    }

    /* {
       #warning remove this
       // test: force fixed vecRMS
       for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
       vecRMS[k] = 0.2; // in arcsec
       }
       }
    */
    
    const double chisq_50  = gsl_cdf_chisq_Pinv(0.50,allOpticalObs.size());
    const double chisq_90  = gsl_cdf_chisq_Pinv(0.90,allOpticalObs.size());
    const double chisq_95  = gsl_cdf_chisq_Pinv(0.95,allOpticalObs.size());
    // const double chisq_99  = gsl_cdf_chisq_Pinv(0.99,allOpticalObs.size());
#warning TEST!
    const double chisq_99  = gsl_cdf_chisq_Pinv(0.99,allOpticalObs.size());
    const double chisq_999 = gsl_cdf_chisq_Pinv(0.999999,allOpticalObs.size());
    // 
    ORSA_DEBUG("chisq 50\%: %.2f  90\%: %.2f  95\%: %.2f  99\%: %.2f  99.9...\%: %.2f",
               chisq_50,chisq_90,chisq_95,chisq_99,chisq_999);
    
    if (1) {
        
        /***** INPUT *****/
        
        const double minAdaptiveRange = orsa::FromUnits(  0.0,orsa::Unit::AU);
        const double maxAdaptiveRange = orsa::FromUnits(100.0,orsa::Unit::AU);
        
        const unsigned int minAdaptiveSize = 2; // 2
        
        const int randomSeed = getpid(); // 717901;
        
        /***** END INPUT *****/
        
        ORSA_DEBUG("randomSeed: %i",randomSeed);
        osg::ref_ptr<orsa::RNG> rng = new orsa::RNG(randomSeed);
        
        if (allOpticalObs.size() < 2) {
            ORSA_DEBUG("problem: not enough observations");
            exit(0);
        }
        
        // observation times must be sorted
        for (unsigned int j=1; j<allOpticalObs.size(); ++j) {
#warning are equal time allowed? that can happen...
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
        
        // std::vector<orsa::Vector> u_o2an; // unit vectors, from Dawn to Satellite, with "noise" (astrometric uncertainty)
        u_o2an.resize(allOpticalObs.size());
        
        // std::vector<orsa::Vector> R_o2an, V_o2an; // real length vectors, from Dawn to Satellite, with "noise" (astrometric uncertainty)
        R_o2an.resize(allOpticalObs.size());
        V_o2an.resize(allOpticalObs.size());
        
        // std::vector<orsa::Vector> R_an, V_an;
        R_an.resize(allOpticalObs.size());
        V_an.resize(allOpticalObs.size());
        
        // std::vector<orsa::Vector> R_s2an;
        R_s2an.resize(allOpticalObs.size());
        
        // std::vector<orsa::Vector> V_s2an;
        V_s2an.resize(allOpticalObs.size());
        
        // computed ONCE ONLY, otherwise it rotates...
        // const orsa::Matrix l2g = orsa::localToGlobal(vesta.get(),bg.get(),obsTime[0]);
        // const orsa::Matrix g2l = orsa::globalToLocal(vesta.get(),bg.get(),obsTime[0]);
        
        for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
            
            const orsa::Time t = allOpticalObs[j]->epoch;
            
            if (!bg->getInterpolatedPosVel(R_s[j],V_s[j],sun.get(),t)) {
                ORSA_DEBUG("problems... t:");
                orsa::print(t);
            }
            
            if (!obsPosCB->getPosVel(R_o[j],V_o[j],allOpticalObs[j].get())) { 
                ORSA_DEBUG("problems... t:");
                orsa::print(t);
            }
            
            u_o2a[j] = observationDirection(allOpticalObs[j]);
            
            xS_o2a[j] = orsa::externalProduct(orsa::Vector(0,0,1),u_o2a[j]).normalized();
            yS_o2a[j] = orsa::externalProduct(u_o2a[j],xS_o2a[j]).normalized();
        } 
        
        {
            // test
            for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                for (unsigned int k=j+1; k<allOpticalObs.size(); ++k) {
                    const double angle = acos(u_o2a[j]*u_o2a[k]);
                    const double dt = (allOpticalObs[k]->epoch-allOpticalObs[j]->epoch).get_d();
                    if (k-j==1) {
                        ORSA_DEBUG("angle between [%02i] and [%02i] obs. = %8.3f [arcsec]   dt = %8.1f [s]   rate = %8.3f [arcsec/s]",
                                   j,
                                   k,
                                   angle*orsa::radToArcsec(),
                                   orsa::FromUnits(dt,orsa::Unit::SECOND,-1),
                                   orsa::FromUnits(angle/dt,orsa::Unit::SECOND)*orsa::radToArcsec());
                    } else {
                        ORSA_DEBUG("angle between [%02i] and [%02i] obs. = %8.3f [arcsec]   dt = %8.1f [s]",
                                   j,
                                   k,
                                   angle*orsa::radToArcsec(),
                                   orsa::FromUnits(dt,orsa::Unit::SECOND,-1));
                    }   
                }
            }
        }
        
        std::vector< orsa::Cache<double> > vecMaxRange;
        vecMaxRange.resize(allOpticalObs.size());
        {
            // test: determine the maximum range for each observation
            osg::ref_ptr<orsaUtil::MultiminMaximumRange> mmr = new orsaUtil::MultiminMaximumRange;
            
            for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
                    if (j==k) continue;
                    const double maxRange = mmr->getMaximumRange(R_s[j],R_o[j],u_o2a[j],vecRMS[j],allOpticalObs[j]->epoch,
                                                                 R_s[k],R_o[k],u_o2a[k],vecRMS[k],allOpticalObs[k]->epoch,
                                                                 orsaSolarSystem::Data::GMSun(),
                                                                 2.5);
#warning nominal max error here... 3.0? 2.5?
                    
                    // maxRange is only for first set of arguments == [j]
                    // vecMaxRange[j].setIfLarger(maxRange);
                    if (maxRange > 0.0) {
                        vecMaxRange[j].setIfLarger(maxRange);
                    }
                    
                    /* ORSA_DEBUG("max range distance for obs [%02i] and [%02i] = %8.3g [AU]",
                       j,
                       k,
                       orsa::FromUnits(maxRange,orsa::Unit::AU,-1));
                    */
                    
                }
                // ORSA_DEBUG("set: %i",vecMaxRange[j].isSet());
                if (vecMaxRange[j].isSet()) {
                    ORSA_DEBUG("max range distance for obs [%02i] = %8.3f [AU]",
                               j,
                               orsa::FromUnits(vecMaxRange[j],orsa::Unit::AU,-1));
                } else {
                    // ORSA_DEBUG("vecMaxRange[%i] not set??",j);
#warning correct to force it here?
                    vecMaxRange[j] = maxAdaptiveRange; 
                }
            }
            
            if (0) {
#warning remove this
                // test: force large maxRange...
                for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                    vecMaxRange[j] = orsa::FromUnits(1.00,orsa::Unit::AU);
                }
            }
            
            // try to reduce vecMaxRange values:
            // take the smallest, and reduce all others according to V_escape*dt....
            orsa::Cache<unsigned int> indexMin;
            // first, find one initial indexMin, any one works, the first that is set
            for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
                if (vecMaxRange[k].isSet()) {
                    indexMin=k;
                    break;
                }
            }
            if (indexMin.isSet()) {
                
                for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
                    if (vecMaxRange[k].isSet()) {
                        if (vecMaxRange[k] < vecMaxRange[indexMin]) indexMin = k;
                    }
                }
                ORSA_DEBUG("SMALLEST max range distance for obs [%02i] = %8.3f [AU]",
                           (*indexMin),
                           orsa::FromUnits(vecMaxRange[indexMin],orsa::Unit::AU,-1));
                
                const orsa::Vector R_a = R_o[indexMin] + u_o2a[indexMin]*vecMaxRange[indexMin];
                const double escapeVelocity = sqrt(2*orsaSolarSystem::Data::GMSun()/(R_a-R_s[indexMin]).length());
                
                for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                    if (j==indexMin) continue;
                    const double rangeCenter = (R_a-R_o[j])*u_o2a[j];
                    
                    const double minDistance = (R_o[j]+u_o2a[j]*rangeCenter-R_a).length();
                    const double dt = fabs((allOpticalObs[j]->epoch-allOpticalObs[indexMin]->epoch).get_d());
                    const double maxDistance = escapeVelocity*dt;
                    if (minDistance > maxDistance) {
                        /* ORSA_DEBUG("------------- minD: %g   maxD: %g",
                           FromUnits(minDistance,orsa::Unit::AU,-1),
                           FromUnits(maxDistance,orsa::Unit::AU,-1));
                        */
                        if (rangeCenter > 0.0) {
                            // now setting if SMALLER!
                            vecMaxRange[j].setIfSmaller(rangeCenter);
                        }
                    } else {
                        const double rangeDelta  = sqrt(orsa::square(maxDistance)-orsa::square(minDistance));
                        const double newRange = rangeCenter+rangeDelta;
                        if (newRange > 0.0) {
                            // now setting if SMALLER!
                            vecMaxRange[j].setIfSmaller(rangeCenter+rangeDelta);
                        }
                    }
                }
                
                for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                    ORSA_DEBUG("REDUCED max range distance for obs [%02i] = %8.3f [AU]",
                               j,
                               orsa::FromUnits(vecMaxRange[j],orsa::Unit::AU,-1));
                }
            }
        }
        
        // keep these in sync!
        // typedef double AdaptiveIntervalTemplateType;
        // typedef Entry AdaptiveIntervalTemplateType;
        // typedef orsaUtil::AdaptiveInterval<AdaptiveIntervalTemplateType> AdaptiveIntervalType;
        // typedef orsaUtil::AdaptiveIntervalElement<AdaptiveIntervalTemplateType> AdaptiveIntervalElementType;
        std::vector< osg::ref_ptr<AdaptiveIntervalType> > range_99;
        range_99.resize(allOpticalObs.size());
        for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
            range_99[k] = new AdaptiveIntervalType(minAdaptiveRange,
                                                   std::min(double(vecMaxRange[k]),maxAdaptiveRange),
                                                   0.99, // confidence level for this interval, different from the chisq-level
                                                   chisq_99,
                                                   randomSeed+234,
                                                   allOpticalObs);
            // #warning RESTORE THIS!? USE vecMaxRange...
            /* range_99[k] = new AdaptiveIntervalType(minAdaptiveRange,
               maxAdaptiveRange,
               0.99,
               chisq_99,
               randomSeed+234);
            */
        }
        
        osg::ref_ptr<orsaUtil::MultiminOrbitalVelocity> mov = new orsaUtil::MultiminOrbitalVelocity;
        unsigned int iter=0;
        // for (unsigned int zzz=0; zzz<10000; ++zzz) {
        //
        // test counts
#warning need to improve this! (one for each value of chisq confidence level?)
        unsigned int ct_tot=0;
        unsigned int ct_NEO=0;
        unsigned int old_ct_tot=0;
        //
#warning remove this initial value of 99999;
        orsa::Cache<double> minRMS = 99999;
        //
        // std::list<Entry> entryList;
        //
        while (iter<100000000) {
            ++iter;
            // ORSA_DEBUG("ITER: %i",iter);
            // first add noise to unit vectors

            // don't use ___z1 and ___z2 in code below, but z1 and z2 
            unsigned int ___z1 = rng->gsl_rng_uniform_int(allOpticalObs.size());
            unsigned int ___z2;
            do { ___z2 = rng->gsl_rng_uniform_int(allOpticalObs.size()); }
            while (___z1 == ___z2);
            // save them as z1 < z2 to ensure that vectors...[z1] are set before vectors...[z2]
            const unsigned int z1 = std::min(___z1,___z2);
            const unsigned int z2 = std::max(___z1,___z2);
            {
                bool bound=true;
                double dx,dy;
                for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                    if ((j!=z1) && (j!=z2)) continue;
                    if ( (j==z1) ||
                         (j==z2 && range_99[z2]->size()>=minAdaptiveSize) ) {
                        // if (j==z2) ORSA_DEBUG("speed-up...");
                        // normal sampling on z1
                        rng->gsl_ran_dir_2d(&dx,&dy);
                        // u_o2an[j] = (u_o2a[j]+rng->gsl_ran_gaussian(astrometricSigma)*(dx*xS_o2a[j]+dy*yS_o2a[j])).normalized();
                        u_o2an[j] = (u_o2a[j]+rng->gsl_ran_gaussian(vecRMS[j]*orsa::arcsecToRad())*(dx*xS_o2a[j]+dy*yS_o2a[j])).normalized();
                        // R_o2an[j] = range->sample()*u_o2an[j];
#warning which range to use?
                        R_o2an[j] = range_99[j]->sample()*u_o2an[j];
                        R_an[j] = R_o[j] + R_o2an[j];
                        R_s2an[j] = R_an[j] - R_s[j];
                        if (j==z2 && range_99[z2]->size()>=minAdaptiveSize) {
                            const double escapeVelocity = sqrt(2*orsaSolarSystem::Data::GMSun()/R_s2an[z1].length());
                            const double dt = (allOpticalObs[z2]->epoch-allOpticalObs[z1]->epoch).get_d();
                            const double maxDistance = escapeVelocity*dt;
                            if ((R_s2an[z2]-R_s2an[z1]).length() > maxDistance) {
                                // no possible bound solutions here
                                bound=false;
                                continue;
                            }
                        }
                    }
                    if (j==z2 && range_99[z2]->size()<minAdaptiveSize) {
                        // enforce bound orbit on z2
                        const double escapeVelocity = sqrt(2*orsaSolarSystem::Data::GMSun()/R_s2an[z1].length());
                        const double dt = (allOpticalObs[z2]->epoch-allOpticalObs[z1]->epoch).get_d();
                        const double maxDistance = escapeVelocity*dt;
                        //
                        rng->gsl_ran_dir_2d(&dx,&dy);
                        // u_o2an[j] = (u_o2a[j]+rng->gsl_ran_gaussian(astrometricSigma)*(dx*xS_o2a[j]+dy*yS_o2a[j])).normalized();
                        u_o2an[j] = (u_o2a[j]+rng->gsl_ran_gaussian(vecRMS[j]*orsa::arcsecToRad())*(dx*xS_o2a[j]+dy*yS_o2a[j])).normalized();
                        //
                        // const orsa::Vector diffVector = R_s2an[z1]-R_o[z2]+R_s[z2];
                        const orsa::Vector diffVector = R_s2an[z1]+R_s[z1]-R_o[z2];
                        const double rangeCenter = std::max(0.0,u_o2an[j]*diffVector); // range of min distance to R...z1
                        const double minDistance = (u_o2an[j]*rangeCenter-diffVector).length();
                        if (minDistance > maxDistance) {
                            // no possible bound solutions here
                            bound=false;
                            continue;
                        }
                        const double rangeDelta  = sqrt(orsa::square(maxDistance)-orsa::square(minDistance));
                        //
                        // ORSA_DEBUG("R_z1: %g  R_z2_center: %g  diff: %g",R_o2an[z1].length(),rangeCenter,R_o2an[z1].length()-rangeCenter);
                        // ORSA_DEBUG("maxDistance: %g   minDistance: %g   rangeDelta: %g",maxDistance,minDistance,rangeDelta);
                        //
                        R_o2an[j] = std::max(0.0,rangeCenter+rangeDelta*(2*rng->gsl_rng_uniform()-1))*u_o2an[j];
                        R_an[j] = R_o[j] + R_o2an[j];
                        R_s2an[j] = R_an[j] - R_s[j];
                    }
                    if (j==z2) break; // done
                }
                if (!bound) continue;
            }
            
            V_s2an[z1] = mov->getOrbitalVelocity(R_s2an[z1],allOpticalObs[z1]->epoch,
                                                 R_s2an[z2],allOpticalObs[z2]->epoch,
                                                 orsaSolarSystem::Data::GMSun());

            {
                orsaSolarSystem::OrbitWithEpoch tmpOrbit;
                tmpOrbit.compute(R_s2an[z1],
                                 V_s2an[z1],
                                 orsaSolarSystem::Data::GMSun());
                tmpOrbit.epoch = allOpticalObs[z1]->epoch;
                const orsaSolarSystem::OrbitWithEpoch O_s2an_g = tmpOrbit;
                
                AdaptiveIntervalElementType e;
                e.O_s2an_g = O_s2an_g;
                 // e.tested   = false;
                // entryList.push_back(e);
#warning CORRECT?
                for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
                    e.position = R_o2an[k].length();
                    range_99[k]->insert(e);
                }
            }
            
            
            // orsa::print(O_v2sn_g);
            
            /* ORSA_DEBUG("FINAL ORBIT: a: %16.6f [km]   e: %8.6f",
               orsa::FromUnits(O_v2sn_g.a,orsa::Unit::KM,-1),
               O_v2sn_g.e);
            */
            
#warning REWRITE THIS PART...
#if 0
            
            std::list<Entry>::const_iterator it = entryList.begin();
            while (it != entryList.end()) {
                
                // perform this check only if enought points are already available, to speed up things...
#warning should replace this threshold (16) with a check that enough points are available with lower RMS value
                if (range_99[0]->size() >= 16) {
                    bool newMinRMS=false;
                    if (!minRMS.isSet()) {
                        minRMS = (*it).RMS;
                        newMinRMS=true;
                    } else {
                        if ((*it).RMS < minRMS) {
                            minRMS = (*it).RMS;
                            newMinRMS=true;
                        }
                    }


                    
                    if (newMinRMS) {
                        ORSA_DEBUG("RESET HERE...  new minRMS = %g   range_99[0]->size(): %i",(*minRMS),range_99[0]->size());
#warning remember to reset all other counters around the code...
                        for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
                            range_99[k]->reset();

                            
                            // reset counters...
                            ct_tot=0;
                            ct_NEO=0;
                            old_ct_tot=0;
                            //
                            if (vec_residual[k]>vecRMSnominal[k]) {
#warning MAXIMUM RMS value set here, with min(...,3.0)...
                                vecRMS[k] = std::min(vec_residual[k],2.5);
                                ORSA_DEBUG("RMS[%02i] = %5.2f (updated)",k,vecRMS[k]);
                            } else {
                                ORSA_DEBUG("RMS[%02i] = %5.2f",k,vecRMS[k]);
                            }
                        }
                        // continue;
                        break;
                    }
                }
                
                {
                    // feedback
                    // #warning create a range for each observation, and cycle on them with the for loop
                    // #warning cannot insert more than 1 in same range! it will bias too much, so need a range for each data point
                    for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                        AdaptiveIntervalElementType e;
                        e.position = R_o2an[j].length();
                        // e.position = R_o2an[0].length();  // use [j] !!!
                        e.level = (*it).chisq;
                        range_99[j]->insert(e);
                    }
                    /* for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                       AdaptiveIntervalElementType e;
                       e.position = V_o2an[j].length();
                       // e.position = R_o2an[0].length();  // use [j] !!!
                       e.level = chisq;
                       rate_99[j]->insert(e);
                       }
                    */
                    // if (iter%1000==0) range_99->print();
                }
                
                {
                    // debug
                    static unsigned int old_range99_0_size=0;
                    if (range_99[0]->size()!=old_range99_0_size) {
                        ORSA_DEBUG("range_99[0]->size(): %i",range_99[0]->size());
                        old_range99_0_size = range_99[0]->size();
                    }
                }
                
                // if (stat_residual->RMS() < successThreshold) {
                if ((*it).chisq < chisq_99) {
                    // stats
                    ++ct_tot;
                    if (O_s2an_g.a*(1.0-O_s2an_g.e) < FromUnits(1.3,orsa::Unit::AU)) {
                        ++ct_NEO;
                    }
                }
                
                (*it).tested=true;
                
                ++it;
            }

#endif // 0 (rewrite part...)
            
#warning test only first?
            if (range_99[0]->size() >= 1000) break;
        }

        ORSA_DEBUG("iter: %i",iter);
        
        // add body to bg
        // #warning fix threshold
        // if (stat_residual->RMS() < successThreshold) {
        /*         if (chisq < chisq_99) {
                   osg::ref_ptr<Body> b = new Body;
                   char bodyName[1024];
                   sprintf(bodyName,"b%i",iter); // fix...
                   b->setName(bodyName);
                   orsa::Vector rOrbit, vOrbit;
                   O_s2an_g.relativePosVel(rOrbit,vOrbit);
                   IBPS ibps;    
                   ibps.time = O_s2an_g.epoch;
                   //
                   ibps.inertial = new PointLikeConstantInertialBodyProperty(0);
                   //
                   ibps.translational = new DynamicTranslationalBodyProperty;
                   #warning make sure epoch(z1) == orbit epoch!
                   if (allOpticalObs[z1]->epoch != ibps.time) {
                   ORSA_DEBUG("problems!");
                   }
                   ibps.translational->setPosition(rOrbit+R_s[z1]);
                   ibps.translational->setVelocity(vOrbit+V_s[z1]);
                   // orsa::print(ibps.translational->position());
                   // orsa::print(ibps.translational->velocity());
                   b->setInitialConditions(ibps);
                   
                   bg->addBody(b.get());
                   
                   ORSA_DEBUG("bg->size(): %i",bg->size());
                   }
        */
        
        // test: slowly decrement all...
        /* {
           range->feedbackAll(0.9999);
           }
        */
            
        // ORSA_DEBUG("RMS: %8.3f",stat_residual->RMS());
        
        AdaptiveIntervalType::DataType::DataType::const_iterator it = range_99[0]->getData()->getData().begin();
        while (it != range_99[0]->getData()->getData().end()) {
            
            // main output
            if ((*it).level < chisq_99) {

                const orsaSolarSystem::OrbitWithEpoch O_s2an_g = (*it).O_s2an_g;
                
                // compute astrometric offset
#warning this code is repeated many times, should write a function...
                std::vector<double> vec_residual;
                vec_residual.resize(allOpticalObs.size());
                orsa::Vector rOrbit, vOrbit;
                for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                    const orsa::Time t = allOpticalObs[j]->epoch;
                    orsaSolarSystem::OrbitWithEpoch tmpOrbit = O_s2an_g;
                    tmpOrbit.M = fmod(O_s2an_g.M + twopi()*(t-O_s2an_g.epoch).get_d()/O_s2an_g.period(),orsa::twopi());
                    tmpOrbit.relativePosVel(rOrbit,vOrbit);
                    rOrbit += R_s[j];
                    vOrbit += V_s[j];
                    // debug: print R_sn[j] before changing it
                    /* ORSA_DEBUG("distance: %g [km]",
                       orsa::FromUnits((R_sn[j]-rOrbit).length(),orsa::Unit::KM,-1));
                    */
                    // note how some vectors are getting redefined
                    R_an[j] = rOrbit;
                    V_an[j] = vOrbit;
                    R_o2an[j] = (R_an[j]-R_o[j]);
                    V_o2an[j] = (V_an[j]-V_o[j]);
                    u_o2an[j] = (R_an[j]-R_o[j]).normalized();
                    const double residual = acos(std::min(1.0,u_o2an[j]*u_o2a[j]))*orsa::radToArcsec();
                    vec_residual[j] = residual;
                }
                
                // a string to write all the single residuals
                char line_eachResidual[1024*1024];
                {
                    char tmpStr[1024];
                    line_eachResidual[0] = '\0';
                    for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                        sprintf(tmpStr,"%6.2f ",vec_residual[j]);
                        strcat(line_eachResidual,tmpStr);
                    }
                }
                
                // this is relative to the observer
                char line_eachDistance[1024*1024];
                {
                    char tmpStr[1024];
                    line_eachDistance[0] = '\0';
                    for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                        sprintf(tmpStr,"%6.3f ",orsa::FromUnits(R_o2an[j].length(),orsa::Unit::AU,-1));
                        strcat(line_eachDistance,tmpStr);
                    }
                }
                
                // this is relative to the Sun
                char line_eachVelocity[1024*1024];
                {
                    char tmpStr[1024];
                    line_eachVelocity[0] = '\0';
                    for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                        sprintf(tmpStr,"%4.1f ",orsa::FromUnits(orsa::FromUnits((V_an[j]-V_s[j]).length(),orsa::Unit::KM,-1),orsa::Unit::SECOND));
                        strcat(line_eachVelocity,tmpStr);
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
                        if (allOpticalObs[j]->mag.isSet()) {
#warning see the correction here to the reported magnitude, depending on the filter 
                            stat_H->insert(absoluteMagnitude(allOpticalObs[j]->mag+orsaSolarSystem::MPC_band_correction(allOpticalObs[j]->band),
                                                             0.15,
                                                             acos((R_o[j]-R_an[j]).normalized()*(R_s[j]-R_an[j]).normalized()),
                                                             R_o2an[j].length(),
                                                             (R_an[j]-R_s[j]).length()));
                        }
                    }
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
                               randomSeed+88,
                               16,
                               1.0e-6);
                }
                
                
                ORSA_DEBUG("SAMPLE[minRMS=%.3f]: %6.3f %7.3f %9.3f %8.6f %7.3f %5.2f %6.4f %s %s %s %s",
                           (*minRMS),
                           (*(*it).RMS),
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
