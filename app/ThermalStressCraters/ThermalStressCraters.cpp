#include "ThermalStressCraters.h"

int main (int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    /*
       {
       // test
       const double heliocentricDistance = orsa::FromUnits(1.0,orsa::Unit::AU);
       double solarCenterElevation = -1.0*orsa::degToRad();
       while (solarCenterElevation<1.0*orsa::degToRad()) {
       const double solarFraction = SolarDiskFraction( heliocentricDistance, solarCenterElevation);
       ORSA_DEBUG("solar elevation: %.9g [arcsec]   solar disk fraction: %.9g",solarCenterElevation*orsa::radToArcsec(),solarFraction);
       solarCenterElevation += 1.0*orsa::arcsecToRad();
        }
        }
    */
    
    /* if (argc != 9999) {
       ORSA_DEBUG("Usage: %s <D,km> <d,km> <lat,deg> <lon,deg> <R,km> <slope-out,deg> <slope-azimuth,deg> <slope-rim,deg> <shape-par> <pole-R.A.,deg> <pole-Dec.,deg> <pole-phi-J2000,deg> <a,AU> <ecc> <i,deg> <node,deg> <peri,deg> <off-North,km> <off-East,km>",argv[0]);
       // exit(0);
       }
    */
    //
    if (argc != 3) {
        ORSA_DEBUG("Usage: %s <dX,km> <dY,km>",argv[0]);
        exit(0);
    }
    
    const double crater_pX = orsa::FromUnits(atof(argv[1]),orsa::Unit::KM);
    const double crater_pY = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const double crater_pR = sqrt(orsa::square(crater_pX)+orsa::square(crater_pY));
    const double crater_phi = atan2(crater_pY,crater_pX);
    
    if (0) {
        // test
        double r=0.0;
        double h,dhdr;
        while (r<=50.0) {
            CraterShape(h,dhdr,r,100,10,tan(0.0*orsa::degToRad()),tan(40.0*orsa::degToRad()));
            ORSA_DEBUG("cs %g %g %g",r,h,atan(dhdr)*orsa::radToDeg());
            r+=0.01;
        }
    }
    

    // body orientation, g2l of sun direction, scalar product with vector normal to central crater point (or to local zenith/radial direction?), if positive then start to check if given point is illuminated...
    
    
    // sun must be above "global" horizon (spherical body) and "local" horizon (crater local tilted reference plane)
    
    
    // note: body position + radius to crater =~ body position alone, when the sun is so far away...
    
    
    // body orbit
    orsa::Orbit orbit;
    orbit.mu = orsaSolarSystem::Data::GMSun();
    orbit.a = orsa::FromUnits(2.54,orsa::Unit::AU);
    orbit.e = 0.00;
    orbit.i = 5.00*orsa::degToRad();
    orbit.omega_node       = 0.0*orsa::degToRad();
    orbit.omega_pericenter = 0.0*orsa::degToRad(); 
    // all values of orbit.M are sampled, to compute Fs
    const double orbit_period = orbit.period();
    
    const double craterDiameter = orsa::FromUnits(20.0,orsa::Unit::KM);
    const double craterDepth    = orsa::FromUnits( 3.0,orsa::Unit::KM);
    const double craterCenterSlope = tan( 0.0*orsa::degToRad());
    const double craterRimSlope    = tan(30.0*orsa::degToRad());
    const double craterLatitude = 45.0*orsa::degToRad();
    // longitude is not relevant, assuming 0.0;
    const double bodyRadius     = orsa::FromUnits(100.0,orsa::Unit::KM);
#warning check that body radius >> crater diameter...
    const double craterPlaneSlope = tan(15.0*orsa::degToRad());
    const double craterPlaneSlopeAzimuth   = 220.0*orsa::degToRad(); // 0=N, 90=E, 180=S, 270=W, points from high to low
    const double bodyPoleEclipticLatitude  =  60.0*orsa::degToRad();
    const double bodyPoleEclipticLongitude =  20.0*orsa::degToRad();
    
    // #warning OFF-North (km) and OFF-East (km) should be arguments!
    
#warning print obliquity...
    
#warning in most cases, the value of bodyRadius is not essential, only the orientation of the crater matters
    
    ORSA_DEBUG("thermal inertia: %g (SI)",thermalInertia());
    
#warning need to double-check these vectors when working at latitudes below equator
    
#warning the crater is surrounded by a plane... should it be on a sphere? more complex, real gain?
#warning also because in general the body radius seems to be not central to the main computations
    
    const orsa::Vector local_u_pole(0,0,1);
    const orsa::Vector local_u_radial(cos(craterLatitude),0.0,sin(craterLatitude));
    const orsa::Vector local_u_east(orsa::externalProduct(local_u_pole,local_u_radial).normalized());
    const orsa::Vector local_u_north(orsa::externalProduct(local_u_radial,local_u_east).normalized());

    // local_u_horizontal is a vector on the h=z=0 level, horizontal (no slope, no matter what the slope of the crater is)
    const orsa::Vector local_u_horizontal =
        cos(craterPlaneSlopeAzimuth-orsa::halfpi())*local_u_north +
        sin(craterPlaneSlopeAzimuth-orsa::halfpi())*local_u_east;
    // local_u_up is up from the crater (zenith?)
    const orsa::Vector local_u_up = orsa::Matrix::axisRotation(local_u_horizontal,craterPlaneSlope)*local_u_radial;
    const orsa::Vector local_u_low_to_high(orsa::externalProduct(local_u_up,local_u_horizontal).normalized());
    
    orsa::print(local_u_radial);
    orsa::print(local_u_east);
    orsa::print(local_u_north);
    orsa::print(local_u_horizontal);
    orsa::print(local_u_up);
    orsa::print(local_u_low_to_high);

    orsa::print(local_u_horizontal*local_u_low_to_high);
    orsa::print(local_u_up*local_u_low_to_high);
    orsa::print(local_u_horizontal*local_u_up);
    
    // so now the 3 cartesial vectos wich are in crater coordinates are local_u_horizontal, local_u_low_to_high (both at h=0), and local_u_up
    
    const orsa::Vector local_crater_center_h0 = local_u_radial*bodyRadius;
    const orsa::Vector local_crater_center    = local_crater_center_h0 - local_u_up*craterDepth;
    //
    double h,dhdr;
    const bool goodCrater = CraterShape(h,
                                        dhdr,
                                        crater_pR,
                                        craterDiameter,
                                        craterDepth,
                                        craterCenterSlope,
                                        craterRimSlope);
    ORSA_DEBUG("cs %g %g %g",crater_pR,h,atan(dhdr)*orsa::radToDeg());
    if (!goodCrater) {
        ORSA_DEBUG("problems...");
        exit(0);
    }
    const double crater_h    = h; // negative in crater
    const double crater_dhdr = dhdr;
    //
    const double crater_point_slope = dhdr;
    const double crater_point_slope_angle = atan(dhdr);
    /* const orsa::Vector local_crater_point_h0 =
       local_crater_center_h0 +
       crater_pX*( cos(crater_phi)*local_u_horizontal + sin(crater_phi)*local_u_low_to_high) +
       crater_pY*(-sin(crater_phi)*local_u_horizontal + cos(crater_phi)*local_u_low_to_high);   const orsa::Vector local_crater_point_h0 =
    */
    const orsa::Vector local_crater_point_h0 =
        local_crater_center_h0 +
        crater_pR*(cos(crater_phi)*local_u_horizontal + sin(crater_phi)*local_u_low_to_high);
    const orsa::Vector local_crater_point = local_crater_point_h0 + local_u_up*h;
    const orsa::Vector local_crater_point_normal =
        cos(crater_point_slope_angle)*local_u_up - sin(crater_point_slope_angle)*(cos(crater_phi)*local_u_horizontal + sin(crater_phi)*local_u_low_to_high);
#warning above: u_radial or u_up ?
    
    orsa::print(local_crater_point);
    orsa::print(local_crater_point_normal);
    
    const double local_crater_point_lat = asin(local_crater_point.getZ()/local_crater_point.length());
    const double local_crater_point_lon = atan2(local_crater_point.getY(),local_crater_point.getX());
    
    const orsa::Time t0_Time(0);
    const double t0 = t0_Time.get_d();
    
    osg::ref_ptr<orsaSolarSystem::ConstantZRotationEcliptic_RotationalBodyProperty> rot =
        new orsaSolarSystem::ConstantZRotationEcliptic_RotationalBodyProperty(t0_Time,
                                                                              0.0,
                                                                              omega(),
                                                                              bodyPoleEclipticLongitude,
                                                                              bodyPoleEclipticLatitude);
    
    History history;
    
#warning should be careful here between sidereal and solar rotation period
    
    /* const size_t NS=10000000;
       const size_t days=ceil(orbit_period/rotationPeriod());
       const double total_simulation_time = days*rotationPeriod(); // slightly larger than orbit_period because of the ceil(...) above (need integer days)
       const double dt = total_simulation_time/NS;
    */
    //
    const orsa::Time dt_Time(mpz_class(1000000)*3);
    orsa::print(dt_Time);
    const double dt = dt_Time.get_d();
    const double solarRotationPeriod = rotationPeriod()/(1.0-rotationPeriod()/orbit_period);
    const size_t rotations = 20; // ceil(orbit_period/solarRotationPeriod); // ceil(orbit_period/rotationPeriod());
    const double total_simulation_time = rotations*solarRotationPeriod;
    const size_t NS = total_simulation_time/dt;
    //
    ORSA_DEBUG("sidereal rotation period: %g [h]   solar rotation period: %g [h]   orbital rotation period: %g [year]",
               orsa::FromUnits(rotationPeriod(),orsa::Unit::HOUR,-1),
               orsa::FromUnits(solarRotationPeriod,orsa::Unit::HOUR,-1),
               orsa::FromUnits(orbit_period,orsa::Unit::YEAR,-1));
    // const double hdist = orsa::FromUnits(3.0,orsa::Unit::AU); // heliocentric distance
    // ORSA_DEBUG("big-theta: %g",theta(hdist));
    std::vector<double> Fs;
    Fs.resize(NS);
    orsa::Vector orbitPosition;
    for (size_t p=0; p<Fs.size(); ++p) {
        const orsa::Time t_Time = dt_Time*p;
        const double t = t_Time.get_d();
        orbit.M = orsa::twopi()*(t_Time-t0_Time).get_d()/orbit_period;
        orbit.relativePosition(orbitPosition);
        const double hdist = orbitPosition.length();
        // ORSA_DEBUG("hdist: %g",orsa::FromUnits(hdist,orsa::Unit::AU,-1));

        
        /* Fs[p] = solar()/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2) *
           std::max(cos(days*orsa::twopi()*(double)p/(double)NS),0.0);
        */
        
        rot->update(t_Time);
        const orsa::Matrix m = orsa::QuaternionToMatrix(rot->getQ());
        orsa::Matrix m_tr;
        orsa::Matrix::transpose(m, m_tr);
        const orsa::Matrix globalToLocal = m_tr;
        // orsa::print(globalToLocal);
        const orsa::Vector local_u_sun = globalToLocal*(-orbitPosition.normalized());
        if (local_u_sun*local_u_up > 0.0) {
            
            //  bool illuminated=true;
            double solarDiskFraction=1.0;
            // this 'if' is needed for ponts outside crater...
            if (crater_pR<=0.5*craterDiameter) {
                // this checks for self-shadowing    
                // const double d = (local_crater_center_h0-local_crater_point)*local_u_up;
                // ORSA_DEBUG("%g %g",crater_h,d);
                const double beta = acos(local_u_up*local_u_sun);
                // const double l = d/cos(beta);
                const double l = fabs(crater_h/cos(beta));
                const orsa::Vector P = local_crater_point+local_u_sun*l; // intersection at h=0
                // solve for point on rim (C1 or C2) in the same direction as P,
                // to estimate solar elevation and possible disk fraction obstructed by crater rim
                // work on h0 plane
                const double a = crater_pR;
                const double b = 0.5*craterDiameter;
                const double cos_delta =
                    (local_crater_center_h0-local_crater_point_h0).normalized() *
                    (P-local_crater_point_h0).normalized();
                const double root = sqrt(a*a*cos_delta*cos_delta-a*a+b*b);
                // ORSA_DEBUG("root: %g",root);
                // here the solution is always the larges one, the only one with x positive
                // const double x1 = a*cos_delta-root;
                const double x2 = a*cos_delta+root;
                // ORSA_DEBUG("x1: %g   x2: %g",x1,x2);
                // orsa::Vector C1 = local_crater_point_h0 + x1*(P-local_crater_point_h0).normalized();
                orsa::Vector C2 = local_crater_point_h0 + x2*(P-local_crater_point_h0).normalized();
                // ORSA_DEBUG("C1: %g",(C1-local_crater_center_h0).length());
                // ORSA_DEBUG("C2: %g",(C2-local_crater_center_h0).length());
                //
                const double solarCenterElevation =
                    acos(local_u_up*(C2-local_crater_point).normalized()) - acos(local_u_up*local_u_sun);
                
                solarDiskFraction = SolarDiskFraction(hdist,solarCenterElevation);
                
                if ( (solarDiskFraction < 0.0) || 
                     (solarDiskFraction > 1.0) ) {
                    ORSA_DEBUG("problems...");
                    exit(0);
                }
                
                // ORSA_DEBUG("(P-local_crater_center_h0).length(): %g   solarCenterElevation: %g [deg]",(P-local_crater_center_h0).length(),solarCenterElevation*orsa::radToDeg());
                
                /* const double dist = (P-local_crater_center_h0).length();
                   if (dist>0.5*craterDiameter) {
                   illuminated=false;
                   }
                */
            }
            
#warning need to compute solar disk fraction also for points on plane outside crater...
            
            // ORSA_DEBUG("d: %g   beta: %g   l: %g   dist: %g",d,beta,l,dist);

            if (solarDiskFraction>0.0) {
                Fs[p] = solar()/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2) *
                    solarDiskFraction *
                    std::max(0.0,local_u_sun*local_crater_point_normal);
            } else {
                Fs[p] = 0.0; // self-shadowing of rim over the crater point
            }
            
        } else {
            Fs[p] = 0.0;
        }       
        
            
        // test
        /* Fs[p] = solar()/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2) *
           std::max(cos(days*orsa::twopi()*(double)p/(double)NS),0.0) *
           std::max(0.5+0.7*cos(orsa::twopi()*(double)p/(double)NS),0.0);
           if (p%99==0) Fs[p] = 0.0;
        */
        
        // TEST!
        /* double proj = std::max(cos(days*orsa::twopi()*(double)p/(double)NS),0.0);
           if (proj < 0.5) proj = 0.0;
           Fs[p] = solar()/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2)*proj;
        */
    }
    
    /* for (size_t p=0; p<NS; ++p) {
       ORSA_DEBUG("Fs[%06i] = %g",p,Fs[p]);
       }
    */
    
    ORSA_DEBUG("ls: %g",skinDepth());
    
    const size_t numSlices=100;
    const double dx = 0.05*skinDepth(); // or a fraction of skinDepth = ls
    // const double dt = days*rotationPeriod()/NS;
    const unsigned int history_skip = 1;
    const double stability_eps   = 1.0e-3;
    const double convergence_eps = 1.0e-6;
    ComputePeriodicThermalHistory(history,
                                  numSlices,
                                  200.0,
                                  Fs,
                                  dx,
                                  dt,
                                  history_skip,
                                  stability_eps,
                                  convergence_eps);
    
    char filename[4096];
    sprintf(filename,"TSC_history_%+.3f_%+.3f.out",orsa::FromUnits(crater_pX,orsa::Unit::KM,-1),orsa::FromUnits(crater_pY,orsa::Unit::KM,-1));
    FILE * fp = fopen(filename,"w");
    for (unsigned int j=0; j<history.size(); ++j) {
        
        unsigned int jpp = (j==(history.size()-1)) ? (0) : (j+1);
        const double dTdt = (history[jpp][0].T-history[j][0].T)/(history_skip*dt);
        
#warning for Fs should write the daily duty cycle, and the max. value of Fs daily
        
        gmp_fprintf(fp,
                    "%9.f %7.3f %7.3f %.3f %+9.6f %+7.3f %+8.3f %+8.3f %+8.3f %+8.3f %8.3f %8.3f %6.3f\n",
                    j*history_skip*dt,
                    Fs[j*history_skip],
                    history[j][0].T,
                    history[j][history[j].size()-2].T, // -2 because -1 is never changed by thermal algo
                    dTdt,
                    orsa::FromUnits(dTdt,orsa::Unit::MINUTE),
                    orsa::FromUnits(crater_pX,orsa::Unit::KM,-1),
                    orsa::FromUnits(crater_pY,orsa::Unit::KM,-1),
                    local_crater_point_lon*orsa::radToDeg(),
                    local_crater_point_lat*orsa::radToDeg(),
                    orsa::FromUnits(local_crater_point.length(),orsa::Unit::KM,-1),
                    orsa::FromUnits(crater_h,orsa::Unit::KM,-1),
                    crater_point_slope_angle*orsa::radToDeg());
        
    }
    fclose(fp);
    
    const double threshold_dTdt = orsa::FromUnits(2,orsa::Unit::MINUTE,-1);
    // orsa::print(threshold_dTdt);
    
    double max_dTdt=0.0, min_dTdt=0.0, max_depth_above_threshold_dTdt=0.0;
    for (unsigned int j=0; j<history.size(); ++j) {

        for (unsigned int k=0; k<history[0].size(); ++k) {
            
            unsigned int jpp = (j==(history.size()-1)) ? (0) : (j+1);
            const double dTdt = (history[jpp][k].T-history[j][k].T)/(history_skip*dt);
            
            if (dTdt>max_dTdt) max_dTdt = dTdt;
            if (dTdt<min_dTdt) min_dTdt = dTdt;
            
            if (fabs(dTdt)>threshold_dTdt) {
                if (dx*k>max_depth_above_threshold_dTdt) {
                    max_depth_above_threshold_dTdt = dx*k;
                }
            }
        }
    }
    //
    sprintf(filename,"TSC_resume_%+.3f_%+.3f.out",orsa::FromUnits(crater_pX,orsa::Unit::KM,-1),orsa::FromUnits(crater_pY,orsa::Unit::KM,-1));
    fp = fopen(filename,"w");
    gmp_fprintf(fp,"%+8.3f %+8.3f %+8.3f %+8.3f %8.3f %8.3f %6.3f  %+7.3f %+7.3f %8.6f\n",
                orsa::FromUnits(crater_pX,orsa::Unit::KM,-1),
                orsa::FromUnits(crater_pY,orsa::Unit::KM,-1),
                local_crater_point_lon*orsa::radToDeg(),
                local_crater_point_lat*orsa::radToDeg(),
                orsa::FromUnits(local_crater_point.length(),orsa::Unit::KM,-1),
                orsa::FromUnits(crater_h,orsa::Unit::KM,-1),
                crater_point_slope_angle*orsa::radToDeg(),
                orsa::FromUnits(max_dTdt,orsa::Unit::MINUTE),
                orsa::FromUnits(min_dTdt,orsa::Unit::MINUTE),
                orsa::FromUnits(max_depth_above_threshold_dTdt,orsa::Unit::METER,-1));
    fclose(fp);
    
    return 0;
}

