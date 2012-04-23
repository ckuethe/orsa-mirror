#include "ThermalStressCraters.h"

int main (int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
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
    orbit.a = orsa::FromUnits(3.54,orsa::Unit::AU);
    orbit.e = 0.20;
    orbit.i = 5.00*orsa::degToRad();
    orbit.omega_node       = 0.0*orsa::degToRad();
    orbit.omega_pericenter = 0.0*orsa::degToRad(); 
    // all values of orbit.M are sampled, to compute Fs
    const double orbit_period = orbit.period();
    
    const double craterDiameter = orsa::FromUnits(100.0,orsa::Unit::KM);
    const double craterDepth    = orsa::FromUnits( 10.0,orsa::Unit::KM);
    const double craterCenterSlope = tan( 0.0*orsa::degToRad());
    const double craterRimSlope    = tan(30.0*orsa::degToRad());
    const double craterLatitude = 55.0*orsa::degToRad();
    // longitude is not relevant, assuming 0.0;
    const double bodyRadius     = orsa::FromUnits(500.0,orsa::Unit::KM);
#warning check that body radius >> crater diameter...
    const double craterPlaneSlope = tan(15.0*orsa::degToRad());
    const double craterPlaneSlopeAzimuth   = 220.0*orsa::degToRad(); // 0=N, 90=E, 180=S, 270=W, points from high to low
    const double bodyPoleEclipticLatitude  =  60.0*orsa::degToRad();
    const double bodyPoleEclipticLongitude =   0.0*orsa::degToRad();
    
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
    //
    const double crater_point_slope = dhdr;
    const double crater_point_slope_angle = atan(dhdr);
    const orsa::Vector local_crater_point =
        local_crater_center_h0 + crater_pX*sin(crater_phi)*local_u_horizontal + crater_pY*cos(crater_phi)*local_u_low_to_high + local_u_up*h;
    const orsa::Vector local_crater_point_normal =
        cos(crater_point_slope_angle)*local_u_radial - sin(crater_point_slope_angle)*(sin(crater_phi)*local_u_horizontal+cos(crater_phi)*local_u_low_to_high);
    
    orsa::print(local_crater_point);
    orsa::print(local_crater_point_normal);
    
    const orsa::Time t0_Time(0);
    const double t0 = t0_Time.get_d();
    
    osg::ref_ptr<orsaSolarSystem::ConstantZRotationEcliptic_RotationalBodyProperty> rot =
        new orsaSolarSystem::ConstantZRotationEcliptic_RotationalBodyProperty(t0_Time,
                                                                              0.0,
                                                                              omega(),
                                                                              bodyPoleEclipticLongitude,
                                                                              bodyPoleEclipticLatitude);
    
    History history;
    const size_t NS=10000000;
    const size_t days=ceil(orbit_period/rotationPeriod());
    const double total_simulation_time = days*rotationPeriod(); // slightly larger than orbit_period because of the ceil(...) above (need integer days)
    const double dt = total_simulation_time/NS;
    const orsa::Time dt_Time(mpz_class(1000000)*dt);
    orsa::print(dt_Time);
    // const double hdist = orsa::FromUnits(3.0,orsa::Unit::AU); // heliocentric distance
    // ORSA_DEBUG("big-theta: %g",theta(hdist));
    std::vector<double> Fs;
    Fs.resize(NS);
    orsa::Vector orbitPosition;
    for (size_t p=0; p<NS; ++p) {
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
            
            bool illuminated=true;
            // if needed for ponts outside crater...
            if (crater_pR<=0.5*craterDiameter) {
                // this checks for self-shadowing    
                const double d = (local_crater_center_h0-local_crater_point)*local_u_up;
                const double beta = acos(local_u_up*local_u_sun);
                const double l = d/cos(beta);
                const orsa::Vector P = local_crater_point+local_u_sun*l;
                const double dist = (P-local_crater_center_h0).length();
                if (dist>0.5*craterDiameter) {
                    illuminated=false;
                }
            }
            
            // ORSA_DEBUG("d: %g   beta: %g   l: %g   dist: %g",d,beta,l,dist);

            if (illuminated) {
                Fs[p] = solar()/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2) *
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
    const double dx = 0.5*skinDepth(); // or a fraction of skinDepth = ls
    // const double dt = days*rotationPeriod()/NS;
    const unsigned int history_skip = 10;
    
    ComputePeriodicThermalHistory(history,
                                  numSlices,
                                  200.0,
                                  Fs,
                                  dx,
                                  dt,
                                  history_skip);
    
    
    FILE * fp = fopen("TSC.out","w");
    for (unsigned int j=0; j<history.size(); ++j) {

        /* unsigned int jmm = (j==0) ? (history.size()-1) : (j-1);
           const double dTdt = (history[j][0].T-history[jmm][0].T)/dt;
        */
        
        unsigned int jpp = (j==(history.size()-1)) ? (0) : (j+1);
        const double dTdt = (history[jpp][0].T-history[j][0].T)/dt;
        
#warning for Fs should write the daily duty cycle, and the max. value of Fs daily
        
        gmp_fprintf(fp,
                    "%g %g %g %g %g\n",
                    j*history_skip*dt,
                    Fs[j*history_skip],
                    history[j][0].T,
                    history[j][history[j].size()-2].T, // -2 because -1 is never changed by thermal algo
                    dTdt);
        
    }
    fclose(fp);
    
    return 0;
}
