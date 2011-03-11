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

// #include "splitRange.h"
// #include "movingWindowRange.h"
#include "AdaptiveInterval.h"


using namespace orsa;

// magnitude function
// alpha = solar phase angle = angle Sun-Asteroid-Observer
// G = slope parameter (G ~= 0.15)
double P (const double & alpha, 
          const double & G) {
    // ORSA_DEBUG("P:   alpha = %f",alpha.get_mpf_t());
    const double phi_1 = exp(-3.33*pow(tan(0.5*alpha),0.63));
    const double phi_2 = exp(-1.87*pow(tan(0.5*alpha),1.22));
    /* 
       ORSA_DEBUG("P = %f   alpha: %f   p1: %f   p2: %f",
       -2.5*log10((1.0-G)*phi_1+G*phi_2),
       alpha.get_mpf_t(),
       phi_1,
       phi_2);
    */
    return (-2.5*log10((1.0-G)*phi_1+G*phi_2));
}

double apparentMagnitude(const double & H,
                         const double & G,
                         const double & phaseAngle,
                         const double & neo2obs,
                         const double & neo2sun) {
    
    const double V = H + P(phaseAngle,G) + 
        5*log10(FromUnits(neo2obs,orsa::Unit::AU,-1)*FromUnits(neo2sun,orsa::Unit::AU,-1));
    
    return V;
}


double absoluteMagnitude(const double & V,
                         const double & G,
                         const double & phaseAngle,
                         const double & neo2obs,
                         const double & neo2sun) {
    
    const double H = V - P(phaseAngle,G) - 
        5*log10(FromUnits(neo2obs,orsa::Unit::AU,-1)*FromUnits(neo2sun,orsa::Unit::AU,-1));
    
    return H;
}



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

inline orsa::Vector observationDirection(const orsaSolarSystem::OpticalObservation * obs) {
    double s_ra,  c_ra;
    double s_dec, c_dec;
    orsa::sincos(obs->ra.getRef(),  &s_ra,  &c_ra);
    orsa::sincos(obs->dec.getRef(), &s_dec, &c_dec);
    return orsaSolarSystem::equatorialToEcliptic() *
        orsa::Vector(c_dec*c_ra,
                     c_dec*s_ra,
                     s_dec);
}

// multimin

class MultiminOrbitalVelocity : public orsa::Multimin {
public:
    double fun(const orsa::MultiminParameters * par) const {
        
        if (0) {
            // debug output 
            for (unsigned int k=0; k<par->size(); ++k) {
                ORSA_DEBUG("par[%02i] = [%20s] = %18.8g   step: %18.8g",
                           k,
                           par->name(k).c_str(),
                           par->get(k),
                           par->getStep(k));
            }
        }
        
        orsa::Orbit orbit;
        orbit.compute(R1.getRef(),
                      getVel(par),
                      mu.getRef());
        orbit.M = fmod(orbit.M+dt.getRef()*orsa::twopi()/orbit.period(),orsa::twopi());
        orsa::Vector Rx;
        if (!orbit.relativePosition(Rx)) {
            ORSA_DEBUG("problems...");
        }
        /* ORSA_DEBUG("a: %16.6f [AU]   e: %8.6f   dR: %16.6f [km]",
           orsa::FromUnits(orbit.a,orsa::Unit::AU,-1),
           orbit.e,
           (Rx-R2.getRef()).length());
        */
        // ORSA_DEBUG("i: %g",orbit.i*orsa::radToDeg());
        // orsa::print(orbit);
        return (Rx-R2.getRef()).lengthSquared(); 
    }
public:
    orsa::Vector getOrbitalVelocity(
        const orsa::Vector & R1_,
        const orsa::Time   & t1_,
        const orsa::Vector & R2_,
        const orsa::Time   & t2_,
        const double       mu_) {
        
        R1 = R1_;
        R2 = R2_;
        mu = mu_;
        dt = (t2_-t1_).get_d();
        u_R = R1_.normalized();
        u_L = orsa::externalProduct(R1_,R2_).normalized();
        if (t2_<t1_) u_L = -u_L.getRef();
        u_N = orsa::externalProduct(u_L.getRef(),u_R.getRef()).normalized();
        
        const double circularVelocity = sqrt(mu_/R1_.length());
        // const double velocityStep     = sqrt(2)*circularVelocity;
        const double velocityStep     = 0.01*circularVelocity;
        
        osg::ref_ptr<orsa::MultiminParameters> par = new orsa::MultiminParameters;
        // parameter: velocity V1 at {R1,t1}, two components only: {Vr,Vn}
        // Vr = V along R1 direction, Vn = V along orthogonal direction
        par->insert("Vr",             0.0,velocityStep);
        par->insert("Vn",circularVelocity,velocityStep);
        //
        setMultiminParameters(par.get());

#warning tune pars here
        if (!run_nmsimplex(128,1.0e-2)) {
            // ORSA_WARNING("the search did not converge.");
        }
        
        osg::ref_ptr<const orsa::MultiminParameters> parFinal = getMultiminParameters();
        //
        return getVel(parFinal.get());
    }
protected:
    // utility
    orsa::Vector getVel(const orsa::MultiminParameters * par) const {
        return (par->get("Vr")*u_R.getRef()+par->get("Vn")*u_N.getRef());
    }
protected:
    orsa::Cache<orsa::Vector> R1, R2;
    orsa::Cache<double> mu;
    orsa::Cache<double> dt;
    orsa::Cache<orsa::Vector> u_R, u_N, u_L; // R=Radial; N=Normal to {L,R}; L=Angular momentum
};

//

int main(int argc, char **argv) {
    
    QApplication app(argc, argv);
    
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
        }
        char line[1024];
        char obsCode[1024];
        double RMS;
        while (fgets(line,1024,fp)) {
            sscanf(line,"%s %lf",obsCode,&RMS);
            obsRMS[obsCode] = RMS;
        }        
        fclose(fp);
    }
    // test
    /* ORSA_DEBUG("RMS[703] = %g",obsRMS["703"]);
       ORSA_DEBUG("RMS[XYZ] = %g",obsRMS["XYZ"]);
    */
    //
    std::vector<double> vecRMS; // same index as observations
    
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
        vecRMS.resize(allOpticalObs.size());
        for (unsigned int k=0; k<allOpticalObs.size(); ++k) {
            double RMS = obsRMS[allOpticalObs[k]->obsCode.getRef()];
            if (RMS==0.0) {
                // if not set (that is, equal to 0.0), use default of 0.7
                ORSA_DEBUG("cannot find nominal accuracy for observatory code [%s], please update file obsRMS.dat",allOpticalObs[k]->obsCode.getRef().c_str());
                RMS=0.7;
            }
            //
            /* allOpticalObs[k]->sigma_ra  = RMS*arcsecToRad();
               allOpticalObs[k]->sigma_dec = RMS*arcsecToRad();
            */
            //
            vecRMS[k] = RMS; // in arcsec
        }
    }
    
    const double chisq_90 = gsl_cdf_chisq_Pinv(0.90,allOpticalObs.size());
    const double chisq_95 = gsl_cdf_chisq_Pinv(0.95,allOpticalObs.size());
    const double chisq_99 = gsl_cdf_chisq_Pinv(0.99,allOpticalObs.size());
    // 
    ORSA_DEBUG("chisq 90\%: %.2f  95\%: %.2f  99\%: %.2f",chisq_90,chisq_95,chisq_99);
    
    if (1) {
        
        /***** INPUT *****/
        
        const int randomSeed = 7175590;
        
        // const double satelliteOrbit_a = orsa::FromUnits(5000.0,orsa::Unit::KM);
        // const double satelliteOrbit_e = 0.05;
        // const double satelliteOrbit_i = 5.0*orsa::degToRad();
        
        // const orsa::Time satelliteOrbit_epoch = orsaSolarSystem::gregorTime(2011,7,1,0,0,0,0);
        
        // observation times (sample times for search field #1)
        /* std::vector<orsa::Time> obsTime;
           obsTime.push_back(orsaSolarSystem::gregorTime(2011,7,12, 7,56,0,0));
           obsTime.push_back(orsaSolarSystem::gregorTime(2011,7,12, 8, 6,0,0));
           obsTime.push_back(orsaSolarSystem::gregorTime(2011,7,12,10, 0,0,0));
           obsTime.push_back(orsaSolarSystem::gregorTime(2011,7,12,10,11,0,0));
           obsTime.push_back(orsaSolarSystem::gregorTime(2011,7,12,18,17,0,0));
           obsTime.push_back(orsaSolarSystem::gregorTime(2011,7,12,18,28,0,0));
           obsTime.push_back(orsaSolarSystem::gregorTime(2011,7,12,20,21,0,0));
           obsTime.push_back(orsaSolarSystem::gregorTime(2011,7,12,20,32,0,0));
        */
        
        // astrometric uncertainty
        // const double astrometricSigma =   1.0*orsa::arcsecToRad(); // for "noise"
        // #warning replace with chisq values... which depend on the confidence level set...
        // const double successThreshold =   1.0;
        // const double outputThreshold  =   3.0;
        
        // distance of object from observer
        // const double R_min = orsa::FromUnits(100.0,orsa::Unit::KM);
        // const double R_max = orsa::FromUnits(  3.0,orsa::Unit::AU);
        
        // double sampledRange; // save it for feedback
        /* osg::ref_ptr<Range> range = new Range(orsa::FromUnits(100.0,orsa::Unit::KM),
           orsa::FromUnits( 3.0,orsa::Unit::AU),
           60,
           randomSeed+33);
        */
        //
        /* osg::ref_ptr<Range> range = new Range(orsa::FromUnits(0.00,orsa::Unit::AU),
           orsa::FromUnits(3.00,orsa::Unit::AU),
           100,
           randomSeed+33);
        */
        //
        /* osg::ref_ptr<MovingWindowRange> range_99 = new  MovingWindowRange(orsa::FromUnits(0.00,orsa::Unit::AU),
           orsa::FromUnits(3.00,orsa::Unit::AU),
           0.50,
           1000,
           chisq_99,
           true,
           randomSeed+44);
        */
        //
        // keep these in sync!
        typedef double AdaptiveIntervalTemplateType;
        typedef AdaptiveInterval<AdaptiveIntervalTemplateType> AdaptiveIntervalType;
        typedef AdaptiveIntervalElement<AdaptiveIntervalTemplateType> AdaptiveIntervalElementType;
        osg::ref_ptr<AdaptiveIntervalType> range_99 =
            new AdaptiveIntervalType(orsa::FromUnits(0.00,orsa::Unit::AU),
                                     orsa::FromUnits(3.00,orsa::Unit::AU),
                                     0.95,
                                     chisq_99,
                                     true,
                                     randomSeed+234);
        
        // const double targetSuccessRate = 0.001;
        // const unsigned int targetSuccessCount = 100000;
        
        /***** END INPUT *****/
        
        ORSA_DEBUG("randomSeed: %i",randomSeed);
        osg::ref_ptr<orsa::RNG> rng = new orsa::RNG(randomSeed);
        
        // orsa::Vector rObj,   vObj,   rVesta,   vVesta,   rSun,   vSun,   rOrbit,   vOrbit;
        // orsa::Vector rSun, rObs, rObj;
        
        if (allOpticalObs.size() < 2) {
            ORSA_DEBUG("problem: not enough observations");
            exit(0);
        }
        
        // observation times must be sorted
        for (unsigned int j=1; j<allOpticalObs.size(); ++j) {
            if (allOpticalObs[j]->epoch.getRef() <= allOpticalObs[j-1]->epoch.getRef()) {
                ORSA_DEBUG("problem: observation times not sorted");
                exit(0);
            }
        }
        
        // s = sun, o = observer, a = asteroid
        std::vector<orsa::Vector> R_s, V_s, R_o; // R_a;
        R_s.resize(allOpticalObs.size());
        V_s.resize(allOpticalObs.size());
        R_o.resize(allOpticalObs.size());
        // R_a.resize(allOpticalObs.size());
        
        // std::vector<orsa::Vector> u_d2v; // unit vectors, from Dawn to Vesta
        // u_d2v.resize(allOpticalObs.size());
        
        // old u_d2s
        std::vector<orsa::Vector> u_o2a; // unit vectors, from Dawn to Satellite
        u_o2a.resize(allOpticalObs.size());
        
        // old xS_d2s
        std::vector<orsa::Vector> xS_o2a, yS_o2a; // unit vectors, orthogonal to u_d2s, to model astrometric accuacy
        xS_o2a.resize(allOpticalObs.size());
        yS_o2a.resize(allOpticalObs.size());
        
        // old u_d2sn
        std::vector<orsa::Vector> u_o2an; // unit vectors, from Dawn to Satellite, with "noise" (astrometric uncertainty)
        u_o2an.resize(allOpticalObs.size());
        
        // old R_d2sn
        std::vector<orsa::Vector> R_o2an; // real length vectors, from Dawn to Satellite, with "noise" (astrometric uncertainty)
        R_o2an.resize(allOpticalObs.size());
        
        // std::vector<orsa::Vector> R_v; // position of Vesta
        // R_v.resize(allOpticalObs.size());
        
        // std::vector<orsa::Vector> R_d; // position of Dawn
        // R_d.resize(allOpticalObs.size());
        
        // std::vector<orsa::Vector> R_s; // position of satellite
        // R_s.resize(allOpticalObs.size());
        
        // std::vector<orsa::Vector> R_sn; // position of satellite with noise
        // R_sn.resize(allOpticalObs.size());
        //
        std::vector<orsa::Vector> R_an, V_an;
        R_an.resize(allOpticalObs.size());
        V_an.resize(allOpticalObs.size());
        
        // std::vector<orsa::Vector> V_sn; // velocity of satellite with noise
        // V_sn.resize(allOpticalObs.size());
        
        // std::vector<orsa::Vector> R_v2sn; // position of satellite with respect to Vesta with noise
        // R_v2sn.resize(allOpticalObs.size());
        
        std::vector<orsa::Vector> R_s2an;
        R_s2an.resize(allOpticalObs.size());
        
        // std::vector<orsa::Vector> V_v2sn; // velocity of satellite with respect to Vesta with noise
        // V_v2sn.resize(allOpticalObs.size());
        //
        std::vector<orsa::Vector> V_s2an;
        V_s2an.resize(allOpticalObs.size());
        
        // computed ONCE ONLY, otherwise it rotates...
        // const orsa::Matrix l2g = orsa::localToGlobal(vesta.get(),bg.get(),obsTime[0]);
        // const orsa::Matrix g2l = orsa::globalToLocal(vesta.get(),bg.get(),obsTime[0]);
        
        // initial velocity of reference orbit, for first output line
        // orsa::Vector V_v2s_0;
        
        // astrometric residuals
        std::vector<double> vec_residual;
        vec_residual.resize(allOpticalObs.size());
        osg::ref_ptr< orsa::Statistic<double> > stat_residual = new orsa::Statistic<double>;
        
        for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
            
            const orsa::Time t = allOpticalObs[j]->epoch.getRef();
            
            if (!bg->getInterpolatedPosVel(R_s[j],V_s[j],sun.get(),t)) {
                ORSA_DEBUG("problems... t:");
                orsa::print(t);
            }
            
            if (!obsPosCB->getPosition(R_o[j],allOpticalObs[j].get())) { 
                ORSA_DEBUG("problems... t:");
                orsa::print(t);
            }
            
            u_o2a[j] = observationDirection(allOpticalObs[j]);
            
            xS_o2a[j] = orsa::externalProduct(orsa::Vector(0,0,1),u_o2a[j]).normalized();
            yS_o2a[j] = orsa::externalProduct(u_o2a[j],xS_o2a[j]).normalized();
        } 

        osg::ref_ptr<MultiminOrbitalVelocity> mov = new MultiminOrbitalVelocity;
        unsigned int iter=-1;
        // for (unsigned int zzz=0; zzz<10000; ++zzz) {
        //
        // test counts
#warning need to improve this! (one for each value of chisq confidence level?)
        unsigned int ct_tot=0;
        unsigned int ct_NEO=0;
        // 
        while (1) {
            ++iter;
            // ORSA_DEBUG("ITER: %i",iter);
            // first add noise to unit vectors

            // don't use ___z1 and ___z2 in code below
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
                    if (j==z1) {
                        // normal sampling on z1
                        rng->gsl_ran_dir_2d(&dx,&dy);
                        // u_o2an[j] = (u_o2a[j]+rng->gsl_ran_gaussian(astrometricSigma)*(dx*xS_o2a[j]+dy*yS_o2a[j])).normalized();
                        u_o2an[j] = (u_o2a[j]+rng->gsl_ran_gaussian(vecRMS[j]*orsa::arcsecToRad())*(dx*xS_o2a[j]+dy*yS_o2a[j])).normalized();
                        // R_o2an[j] = range->sample()*u_o2an[j];
#warning which range to use?
                        R_o2an[j] = range_99->sample()*u_o2an[j];
                        R_an[j] = R_o[j] + R_o2an[j];
                        R_s2an[j] = R_an[j] - R_s[j];
                    }
                    if (j==z2) {
                        // enforce bound orbit on z2
                        const double escapeVelocity = sqrt(2*orsaSolarSystem::Data::GMSun()/R_s2an[z1].length());
                        const double dt = (allOpticalObs[z2]->epoch.getRef()-allOpticalObs[z1]->epoch.getRef()).get_d();
                        const double maxDistance = escapeVelocity*dt;
                        //
                        rng->gsl_ran_dir_2d(&dx,&dy);
                        // u_o2an[j] = (u_o2a[j]+rng->gsl_ran_gaussian(astrometricSigma)*(dx*xS_o2a[j]+dy*yS_o2a[j])).normalized();
                        u_o2an[j] = (u_o2a[j]+rng->gsl_ran_gaussian(vecRMS[j]*orsa::arcsecToRad())*(dx*xS_o2a[j]+dy*yS_o2a[j])).normalized();
                        //
                        const orsa::Vector diffVector = R_s2an[z1]-R_o[z2]+R_s[z2];
                        const double rangeCenter = u_o2an[j]*diffVector; // range of min distance to R...z1
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
                }
                if (!bound) continue;
            }
            
            /* 
               if (firstIter) {
               // force nominal solution
               ORSA_DEBUG("first iter!");
               V_v2sn[0] = V_v2s_0;
               for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
               // "remove noise"....
               u_d2sn[j] = u_d2s[j];
               R_d2sn[j] = R_s[j]-R_d[j];
               R_sn[j]   = R_s[j];
               R_v2sn[j] = R_s[j]-R_v[j];
               }
               } else {
            */
            // solve for a velocity that goes trough the first two points
            // osg::ref_ptr<MultiminOrbitalVelocity> mov = new MultiminOrbitalVelocity;
            /* V_v2sn[0] = mov->getOrbitalVelocity(R_v2sn[0],obsTime[0],
               R_v2sn[1],obsTime[1],
               muVesta);
            */
            // test: using first and last observation
            /* V_v2sn[0] = mov->getOrbitalVelocity(R_v2sn[0],obsTime[0],
               R_v2sn[allOpticalObs.size()-1],obsTime[allOpticalObs.size()-1],
               muVesta);
            */
            // even better: choose two randomly!

            /* unsigned int z1 = rng->gsl_rng_uniform_int(allOpticalObs.size());
               unsigned int z2;
               do { z2 = rng->gsl_rng_uniform_int(allOpticalObs.size()); }
               while (z1 == z2);
            */
            //
            
            V_s2an[z1] = mov->getOrbitalVelocity(R_s2an[z1],allOpticalObs[z1]->epoch.getRef(),
                                                 R_s2an[z2],allOpticalObs[z2]->epoch.getRef(),
                                                 orsaSolarSystem::Data::GMSun());
            // ORSA_DEBUG("%i-%i",z1,z2);
            // }
            
            // check that velocity gives a "bound" orbit
            
            // O_v2sn_g orbit is Vesta-centered and in global ref system
            // all orbits defined at t=obsTime[0];
            orsaSolarSystem::OrbitWithEpoch tmpOrbit;
            tmpOrbit.compute(R_s2an[z1],
                             V_s2an[z1],
                             orsaSolarSystem::Data::GMSun());
            tmpOrbit.epoch = allOpticalObs[z1]->epoch.getRef();
            const orsaSolarSystem::OrbitWithEpoch O_s2an_g = tmpOrbit;
            
            // orsa::print(O_v2sn_g);
            
            /* ORSA_DEBUG("FINAL ORBIT: a: %16.6f [km]   e: %8.6f",
               orsa::FromUnits(O_v2sn_g.a,orsa::Unit::KM,-1),
               O_v2sn_g.e);
            */
            
            // compute astrometric offset
            stat_residual->reset();
            double chisq=0.0;
            orsa::Vector rOrbit, vOrbit;
            for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                const orsa::Time t = allOpticalObs[j]->epoch.getRef();
                tmpOrbit = O_s2an_g;
                tmpOrbit.M = fmod(O_s2an_g.M + twopi()*(t-O_s2an_g.epoch.getRef()).get_d()/O_s2an_g.period(),orsa::twopi());
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
                u_o2an[j] = (R_an[j]-R_o[j]).normalized();
                const double residual = acos(std::min(1.0,u_o2an[j]*u_o2a[j]))*orsa::radToArcsec();
                vec_residual[j] = residual;
                stat_residual->insert(residual);
                chisq+=orsa::square(residual/vecRMS[j]);
                // ORSA_DEBUG("OFFSET[%i]: %5.1f",j,residual);
            }
            
            // O_v2sn_l orbit is in local ref. sys.
            /* tmpOrbit.compute(R_s2an[z1],
               V_s2an[z1],
               orsaSolarSystem::Data::GMSun());
               tmpOrbit.epoch = allOpticalObs[z1]->epoch.getRef();
               const orsaSolarSystem::OrbitWithEpoch O_s2an_l = tmpOrbit;
            */
            
            // feedback x Range class
            // from smaller RMS to largest
            /* for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
               const double sampledRange = R_o2an[j].length();
               if (stat_residual->RMS() < astrometricSigma) {
               range->feedback(sampledRange, 2.000);
               } else if (stat_residual->RMS() < successThreshold) {
               range->feedback(sampledRange, 1.100);
               } else {
               // range->feedback(sampledRange, 1.000);
               }
               }
            */
            //
            /*
               {
               // feedback
               #warning create a range for each observation, and cycle on them with the for loop
               for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
               MovingWindowRange::MovingWindowRangeElement e;
               e.value = R_o2an[j].length();
               // e.value = R_o2an[0].length();  // use [j] !!!
               e.level = chisq;
               range_99->insert(e);
               }
               if (iter%1000==0) range_99->print();
               }
            */
            //
            {
                // feedback
#warning create a range for each observation, and cycle on them with the for loop
#warning cannot insert more than 1 in same range! it will bias too much, so need a range for each data point
                // for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                AdaptiveIntervalElementType e;
                // e.position = R_o2an[j].length();
                e.position = R_o2an[0].length();  // use [j] !!!
                e.level = chisq;
                range_99->insert(e);
                // }
                // if (iter%1000==0) range_99->print();
            }
            
            // if (stat_residual->RMS() < successThreshold) {
            if (chisq < chisq_99) {
                // stats
                ++ct_tot;
                if (O_s2an_g.a*(1.0-O_s2an_g.e) < FromUnits(1.3,orsa::Unit::AU)) {
                    ++ct_NEO;
                }
            }
            
            // add body to bg
            // #warning fix threshold
            // if (stat_residual->RMS() < successThreshold) {
            if (chisq < chisq_99) {
                osg::ref_ptr<Body> b = new Body;
                char bodyName[1024];
                sprintf(bodyName,"b%i",iter); // fix...
                b->setName(bodyName);
                orsa::Vector rOrbit, vOrbit;
                O_s2an_g.relativePosVel(rOrbit,vOrbit);
                IBPS ibps;    
                ibps.time = O_s2an_g.epoch.getRef();
                //
                ibps.inertial = new PointLikeConstantInertialBodyProperty(0);
                //
                ibps.translational = new DynamicTranslationalBodyProperty;
#warning make sure epoch(z1) == orbit epoch!
                if (allOpticalObs[z1]->epoch.getRef() != ibps.time.getRef()) {
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
            
            // test: slowly decrement all...
            /* {
               range->feedbackAll(0.9999);
               }
            */
            
            // ORSA_DEBUG("RMS: %8.3f",stat_residual->RMS());
            
            // main output
            // if (stat_residual->RMS() < outputThreshold) {
            if (chisq < chisq_99) {
                
                // a string to write all the single residuals
                char line_eachResidual[1024];
                {
                    char tmpStr[1024];
                    line_eachResidual[0] = '\0';
                    for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                        sprintf(tmpStr,"%6.2f ",vec_residual[j]);
                        strcat(line_eachResidual,tmpStr);
                    }
                }

                // this is relative to the observer
                char line_eachDistance[1024];
                {
                    char tmpStr[1024];
                    line_eachDistance[0] = '\0';
                    for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                        sprintf(tmpStr,"%6.3f ",orsa::FromUnits(R_o2an[j].length(),orsa::Unit::AU,-1));
                        strcat(line_eachDistance,tmpStr);
                    }
                }
                
                // this is relative to the Sun
                char line_eachVelocity[1024];
                {
                    char tmpStr[1024];
                    line_eachVelocity[0] = '\0';
                    for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                        sprintf(tmpStr,"%4.1f ",orsa::FromUnits(orsa::FromUnits((V_an[j]-V_s[j]).length(),orsa::Unit::KM,-1),orsa::Unit::SECOND));
                        strcat(line_eachVelocity,tmpStr);
                    }
                }
                
#warning improve here
                // absolute magnitude (rough way, is there a better approach?)
                // filter/band correction?
                double H;
                {
                    osg::ref_ptr< orsa::Statistic<double> > stat_H = new orsa::Statistic<double>;
                    for (unsigned int j=0; j<allOpticalObs.size(); ++j) {
                        if (allOpticalObs[j]->mag.isSet()) {
                            stat_H->insert(absoluteMagnitude(allOpticalObs[j]->mag.getRef(),
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
                    EarthOrbit.compute(earth.get(),sun.get(),bg.get(),O_s2an_g.epoch.getRef());
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
                
                
                ORSA_DEBUG("SAMPLE: %6.2f %7.3f %9.3f %8.6f %7.3f %5.2f %6.4f %s %s %s",
                           stat_residual->RMS(),
                           chisq,
                           orsa::FromUnits(O_s2an_g.a,orsa::Unit::AU,-1),
                           O_s2an_g.e,
                           O_s2an_g.i*orsa::radToDeg(),
                           H,
                           orsa::FromUnits(EarthMOID,orsa::Unit::AU,-1),
                           line_eachResidual,
                           line_eachDistance,
                           line_eachVelocity);
                
            } else {
                // ORSA_DEBUG("above threshold: %g",stat_residual->RMS());
            }
            
            // firstIter=false;
            
            /* if (adapt_success >= targetSuccessCount){
               break;
               }
            */
            
            if (0) {
                // debug output
                if (iter%10==0) {
                    if (ct_tot!=0) {
                        ORSA_DEBUG("prob. NEO: %7.3f (%i/%i)",(double)ct_NEO/(double)ct_tot,ct_NEO,ct_tot);
                    }
                }
            }
            
            /* 
               if (0) {
               // debug output
               if (iter%10==0) {
               const Range::DataType & data = range->getData();
               for (unsigned int k=0; k<data.size(); ++k) {
               if (data[k].weight > 1.00) {
               ORSA_DEBUG("%5.3f - %5.3f   w = %5.3f",
               orsa::FromUnits(data[k].min,orsa::Unit::AU,-1),
               orsa::FromUnits(data[k].max,orsa::Unit::AU,-1),
               data[k].weight);
               }
               }
               }
               }
            */
            
            // if (bg->size()>=100) break;
            
        }
    }
    
    ORSA_DEBUG("done.");

    ORSA_DEBUG("bg->size(): %i",bg->size());
    
    osg::ref_ptr<orsa::IntegratorRadau> radau = new orsa::IntegratorRadau;
    radau->integrate(bg.get(),
                     allOpticalObs[0]->epoch.getRef(),
                     allOpticalObs[allOpticalObs.size()-1]->epoch.getRef()+orsa::Time(0,0,0,0,0),
                     orsa::Time(0,0,5,0,0));
    
    {
        // viz
        osg::ref_ptr<orsaOSG::Viz> viz = new orsaOSG::Viz(bg.get(),
                                                          3600.0);
        ViewerQT * viewerWindow = new ViewerQT(0);
        osg::DisplaySettings::instance()->setNumMultiSamples(4);
        osg::Group * rootNode = viz->createRoot();
        viz->_at->setCentralBody(bg->getBody("EARTH"));
        osg::ref_ptr<DepthPartitionNode> dpn = new DepthPartitionNode;
        dpn->addChild(rootNode);
        dpn->setActive(true);
        viewerWindow->setSceneData(dpn.get());
        // depth partion node only supports single window/single threaded at present.
        viewerWindow->setThreadingModel(osgViewer::Viewer::SingleThreaded);
        orsaOSG::FindNamedNodeVisitor fnnv("EARTH");
        rootNode->accept(fnnv);
        if (!fnnv._foundNodes.empty()) {
            // set up the node tracker.
            // osgGA::NodeTrackerManipulator * tm = new osgGA::NodeTrackerManipulator;
            VestaNodeTrackerManipulator * tm = new VestaNodeTrackerManipulator;
            osgGA::NodeTrackerManipulator::TrackerMode   trackerMode =
                osgGA::NodeTrackerManipulator::NODE_CENTER; // NODE_CENTER_AND_ROTATION
            osgGA::NodeTrackerManipulator::RotationMode rotationMode =
                osgGA::NodeTrackerManipulator::TRACKBALL;   // TRACKBALL ELEVATION_AZIM
            tm->setTrackerMode(  trackerMode);
            tm->setRotationMode(rotationMode);
            tm->setTrackNode(fnnv._foundNodes.front().get());
            viewerWindow->setCameraManipulator(tm);
        }
        viewerWindow->show();
        // viewerWindow->showMaximized();
        ORSA_DEBUG("viewer visible: %i",viewerWindow->isVisible());
    }
    
    return app.exec();
}
