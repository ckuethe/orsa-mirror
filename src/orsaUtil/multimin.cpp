#include <orsaUtil/multimin.h>

#include <orsa/orbit.h>

double orsaUtil::MultiminOrbitalVelocity::fun(const orsa::MultiminParameters * par) const {
        
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
    orbit.compute(R1,
                  getVel(par),
                  mu);
    orbit.M = fmod(orbit.M+dt*orsa::twopi()/orbit.period(),orsa::twopi());
    orsa::Vector Rx;
    if (!orbit.relativePosition(Rx)) {
        ORSA_DEBUG("problems...");
    }
    /* ORSA_DEBUG("a: %16.6f [AU]   e: %8.6f   dR: %16.6f [km]",
       orsa::FromUnits(orbit.a,orsa::Unit::AU,-1),
       orbit.e,
       (Rx-R2).length());
    */
    // ORSA_DEBUG("i: %g",orbit.i*orsa::radToDeg());
    // orsa::print(orbit);
    return (Rx-R2).lengthSquared(); 
}

orsa::Vector  orsaUtil::MultiminOrbitalVelocity::getOrbitalVelocity(
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
    if (t2_<t1_) u_L = -u_L;
    u_N = orsa::externalProduct(u_L,u_R).normalized();
        
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
        
    if (!run_nmsimplex(128,1.0e-2)) {
        // ORSA_WARNING("the search did not converge.");
    }
        
    osg::ref_ptr<const orsa::MultiminParameters> parFinal = getMultiminParameters();
    //
    return getVel(parFinal.get());
}

// multimin minimum, maximum range

// auxiliary function
void orsaUtil::MultiminMinMaxRange::fun_plain(double & minDistance,
                                              double & sigmaDistance,
                                              double & maxDistance,
                                              const orsa::MultiminParameters * par) const {
    
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

    const orsa::Vector   R_a_1 = R_o_1 + u_o2a_1*par->get("R_o2a_1");
    const orsa::Vector R_s2a_1 = R_a_1 - R_s_1;
    //
    const double escapeVelocity = sqrt(2*GM_s/R_s2a_1.length());
    const double             dt = fabs((epoch_2-epoch_1).get_d());
    //
    maxDistance = escapeVelocity*dt;
    //
    const orsa::Vector diffVector = R_s2a_1-R_o_2+R_s_2;
    const double rangeCenter = std::max(0.0,u_o2a_2*diffVector);
    //
    minDistance = (u_o2a_2*rangeCenter-diffVector).length();
    // correct minDistance for pointing uncertainty
    //
    // sigmaDistance = totalSigma*rangeCenter; // sin(tS) ~= tS because is small
    sigmaDistance = totalSigma*par->get("R_o2a_1");
    
    /* ORSA_DEBUG("sD: %g  sD: %g  tS: %g  rC: %g  R_o2a_1: %g",
       totalSigma*rangeCenter,
       totalSigma*par->get("R_o2a_1"),
       (*totalSigma),
       rangeCenter,
       par->get("R_o2a_1"));
    */
    
    // should return 0.0 if minDistance < sigmaDistance, but return fabs(mD-sD) in order for the minimization algorithm to find the exact position where mD=sD
    // return fabs(mD-sD);
    //
    /* if (minDistance>sigmaDistance) {
       return (minDistance-sigmaDistance);
       } else {
       return 0.0;
       }
    */
}

// main minimization function
double orsaUtil::MultiminMinMaxRange::fun(const orsa::MultiminParameters * par) const {
    double minDistance, sigmaDistance, maxDistance;
    fun_plain(minDistance,sigmaDistance,maxDistance,par);
    // correct minDistance
    double minDistanceSigma = (minDistance>sigmaDistance) ? (minDistance-sigmaDistance) : 0.0;
    // use fabs(...) in order for the minimization algorithm to find the exact position of the minimum
    return (fabs(minDistanceSigma-maxDistance));
}

void orsaUtil::MultiminMinMaxRange::getMinMaxRange(orsa::Cache<double> & minRange,
                                                   orsa::Cache<double> & maxRange,
                                                   const double & minBoundary,
                                                   const double & maxBoundary,
                                                   const orsa::Vector & R_s_1_,
                                                   const orsa::Vector & R_o_1_,
                                                   const orsa::Vector & u_o2a_1_,
                                                   const double       & sigma_1_arcsec_,
                                                   const orsa::Time   & epoch_1_,
                                                   const orsa::Vector & R_s_2_,
                                                   const orsa::Vector & R_o_2_,
                                                   const orsa::Vector & u_o2a_2_,
                                                   const double       & sigma_2_arcsec_,
                                                   const orsa::Time   & epoch_2_,
                                                   const double       & GM_s_,
                                                   const double       & sigma_factor_) {
    R_s_1   = R_s_1_;
    R_o_1   = R_o_1_;
    u_o2a_1 = u_o2a_1_;
    // sigma_1_arcsec = sigma_1_arcsec_;
    epoch_1 = epoch_1_;
    //
    R_s_2   = R_s_2_;
    R_o_2   = R_o_2_;
    u_o2a_2 = u_o2a_2_;
    // sigma_2_arcsec = sigma_2_arcsec_;
    epoch_2 = epoch_2_;
    //
    GM_s    = GM_s_;
    // sigma_factor = sigma_factor_;
    totalSigma = sigma_factor_*(sigma_1_arcsec_+sigma_2_arcsec_)*orsa::arcsecToRad();
    
    if ( (sigma_1_arcsec_ <= 0.0) ||
         (sigma_2_arcsec_ <= 0.0) ||
         (sigma_factor_   <= 0.0) ) {
        ORSA_DEBUG("problems...");
        exit(0);
    }

    // 
    
    
    
    if (0) {
        // test
        osg::ref_ptr<orsa::MultiminParameters> par = new orsa::MultiminParameters;
        par->insert("R_o2a_1",
                    orsa::FromUnits(100.0,orsa::Unit::AU),
                    orsa::FromUnits(  0.001,orsa::Unit::AU));
        double R_o2a_1 = 0.0;
        while (R_o2a_1 < orsa::FromUnits(100.0,orsa::Unit::AU)) {
            par->set("R_o2a_1",R_o2a_1);
            double f = fun(par.get());
            ORSA_DEBUG("TTT %g %g",
                       orsa::FromUnits(R_o2a_1,orsa::Unit::AU,-1),
                       orsa::FromUnits(      f,orsa::Unit::AU,-1));
            R_o2a_1 += orsa::FromUnits(0.1,orsa::Unit::AU);
        }
        
    }

#warning still needed this?
    /* // initial check
       // if (acos(u_o2a_1*u_o2a_2) < sigma_factor*(sigma_1_arcsec+sigma_2_arcsec)*orsa::arcsecToRad()) {
       if (acos(u_o2a_1*u_o2a_2) < totalSigma) {
       // ORSA_DEBUG("obs. uncertainty is larger than their angular separation, skipping...");
       return 0.0;
       }
    */
    
    osg::ref_ptr<orsa::MultiminParameters> par = new orsa::MultiminParameters;
    // initial values dummy for now...
    par->insert("R_o2a_1",
                orsa::FromUnits(100.0,orsa::Unit::AU),
                orsa::FromUnits(  0.01,orsa::Unit::AU));
    par->setRangeMin("R_o2a_1",0.0);
    setMultiminParameters(par.get());

    {
        // roughly search min and max on a fixed step grid
        std::vector<double> fb;
        //
        double b=minBoundary;
        const double db=(maxBoundary-minBoundary)*0.01;
        double minDistance, sigmaDistance, maxDistance;
        double minDistanceSigma;
        //
        do {
            par->set("R_o2a_1",b);
            fun_plain(minDistance,sigmaDistance,maxDistance,par);
            minDistanceSigma = (minDistance>sigmaDistance) ? (minDistance-sigmaDistance) : 0.0;
            fb.push_back(minDistanceSigma-maxDistance); // no fabs here!
            b += db;
        } while (b <= maxBoundary);
        //
        bool doneMin=false, doneMax=false;
        for (unsigned int k=0; k<fb.size()-1; ++k) {
            /* ORSA_DEBUG("fb[%i]: %g   fb[%i] %g",
               k,fb[k],
               k+1,fb[k+1]);
            */
            // min test
            if (!doneMin) {
                // min search if first few fb values are positive and fb[1] < fb[0], i.e.: positive and decreasing...
                if ( (fb[0] > 0.0) &&
                     (fb[k+1] < 0.0) ) {
                    // search for minRange here
                    // ORSA_DEBUG("minRange search...");
                    par->set("R_o2a_1",minBoundary+db*(k+0.5));
                    par->setRange("R_o2a_1",
                                  minBoundary+db*k,
                                  minBoundary+db*(k+1));
                    setMultiminParameters(par.get());
                    if (run_nmsimplex(1024,1.0e-2)) {
                        osg::ref_ptr<const orsa::MultiminParameters> parFinal = getMultiminParameters();
                        minRange = par->get("R_o2a_1");
                        doneMin=true;
                    } else {
                        ORSA_WARNING("the search did not converge.");
                    }
                }
            }
            if (!doneMax) {
                if ( (fb[k] < 0.0) &&
                     (fb[k+1] > 0.0) ) {
                    // search for maxRange here
                    // ORSA_DEBUG("maxRange search...");
                    par->set("R_o2a_1",minBoundary+db*(k+0.5));
                    par->setRange("R_o2a_1",
                                  minBoundary+db*k,
                                  minBoundary+db*(k+1));
                    setMultiminParameters(par.get());
                    if (run_nmsimplex(1024,1.0e-2)) {
                        osg::ref_ptr<const orsa::MultiminParameters> parFinal = getMultiminParameters();
                        maxRange = par->get("R_o2a_1");
                        doneMax=true;
                    } else {
                        ORSA_WARNING("the search did not converge.");
                    }
                }
            }
            if (doneMin && doneMax) break;
        }
    }
    
    /* if (minRange.isSet()) {
       ORSA_DEBUG("minRange: %g [AU]",orsa::FromUnits(minRange,orsa::Unit::AU,-1));
       }
       if (maxRange.isSet()) {
       ORSA_DEBUG("maxRange: %g [AU]",orsa::FromUnits(maxRange,orsa::Unit::AU,-1));
       }
    */
}
