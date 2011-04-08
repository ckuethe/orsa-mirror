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

// multimin maximum range

double orsaUtil::MultiminMaximumRange::fun(const orsa::MultiminParameters * par) const {
        
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
    const double    maxDistance = escapeVelocity*dt;
    //
    const orsa::Vector u_L = (R_a_1 - R_o_2).normalized();
    const double         L = (R_a_1 - R_o_2).length();


    /* if (fabs(maxDistance/L) > 1.0) {
       return 1.0e20;
       }
    */
    // const double     gamma = asin(maxDistance/L);
    const double gamma = (fabs(maxDistance/L)<1.0 ? asin(maxDistance/L) : orsa::pi());
        
    const double  beta = acos(u_o2a_2*u_L);
    //
    const double dBeta = sigma_factor*(sigma_1_arcsec+sigma_2_arcsec)*orsa::arcsecToRad();
    //
    /* if (dBeta > beta) {
       return 1.0e20;
       }
    */
    //
    // const double   newBeta = beta - dBeta;
        
    /* ORSA_DEBUG("beta: %g [deg]   dBeta: %g [deg]   gamma: %g [deg]   R_o2a_1: %g [AU]   retVal: %g",
       beta*orsa::radToDeg(),
       dBeta*orsa::radToDeg(),
       gamma*orsa::radToDeg(),
       orsa::FromUnits(par->get("R_o2a_1"),orsa::Unit::AU,-1),
       orsa::square((beta-dBeta-gamma)*orsa::radToDeg()));
    */
        
        
    /* 
       const double dbg = fabs(beta-gamma);
       if (dBeta > dbg) {
       return orsa::square(dBeta-dbg);
       } else {
       return 
       }
    */
        
    return orsa::square(beta-dBeta-gamma);
        
    // return orsa::square(beta-dBeta-gamma);
        
    /* 
       const double delta = L*fabs(sin(newBeta));
       // const double delta = L*sin(beta);
           
       return orsa::square(delta-maxDistance);
    */
        
        
        
}

double orsaUtil::MultiminMaximumRange::getMaximumRange(const orsa::Vector & R_s_1_,
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
    sigma_1_arcsec = sigma_1_arcsec_;
    epoch_1 = epoch_1_;
    //
    R_s_2   = R_s_2_;
    R_o_2   = R_o_2_;
    u_o2a_2 = u_o2a_2_;
    sigma_2_arcsec = sigma_2_arcsec_;
    epoch_2 = epoch_2_;
    //
    GM_s    = GM_s_;
    sigma_factor = sigma_factor_;

    if ( (sigma_1_arcsec <= 0.0) ||
         (sigma_2_arcsec <= 0.0) ||
         (sigma_factor   <= 0.0) ) {
        ORSA_DEBUG("problems...");
        exit(0);
    }

    // initial check
    if (acos(u_o2a_1*u_o2a_2) < sigma_factor*(sigma_1_arcsec+sigma_2_arcsec)*orsa::arcsecToRad()) {
        ORSA_DEBUG("obs. uncertainty is larger than their angular separation, skipping...");
        return 0.0;
    }
        
    osg::ref_ptr<orsa::MultiminParameters> par = new orsa::MultiminParameters;
    par->insert("R_o2a_1",
                orsa::FromUnits( 1.00,orsa::Unit::AU),
                orsa::FromUnits( 0.01,orsa::Unit::AU));
    par->setRangeMin("R_o2a_1",0.0);
    //
    setMultiminParameters(par.get());
        
    if (!run_nmsimplex(1024,1.0e-2)) {
        // ORSA_WARNING("the search did not converge.");
    }
        
    osg::ref_ptr<const orsa::MultiminParameters> parFinal = getMultiminParameters();
    //
    return par->get("R_o2a_1");
}
