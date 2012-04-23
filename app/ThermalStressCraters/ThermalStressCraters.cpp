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
        ORSA_DEBUG("Usage: %s <off-North,km> <off-East,km>",argv[0]);
        exit(0);
    }
    
    const double delta_North = orsa::FromUnits(atof(argv[1]),orsa::Unit::KM);
    const double delta_East  = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    
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
    orbit.a = orsa::FromUnits(3.2,orsa::Unit::AU);
    orbit.e = 0.00;
    orbit.i = 0.00*orsa::degToRad();
    orbit.omega_node       = 0.0*orsa::degToRad();
    orbit.omega_pericenter = 0.0*orsa::degToRad(); 
    // all values of orbit.M are sampled, to compute Fs
    
    const double craterDiameter = orsa::FromUnits(50.0,orsa::Unit::KM);
    const double craterDepth    = orsa::FromUnits( 8.0,orsa::Unit::KM);
    const double craterCenterSlope = tan( 0.0*orsa::degToRad());
    const double craterRimSlope    = tan(40.0*orsa::degToRad());
    const double craterLatitude = 30.0*orsa::degToRad();
    // longitude is not relevant, assuming 0.0;
    const double bodyRadius     = orsa::FromUnits(500.0,orsa::Unit::KM);
#warning check that body radius >> crater diameter...
    const double craterPlaneSlope = tan(0.0*orsa::degToRad());
    const double craterPlaneSlopeAzimuth   =  0.0*orsa::degToRad(); // 0=N, 90=E, 180=S, 270=W, points from high to low
    const double bodyPoleEclipticLatitude  =  90.0*orsa::degToRad();
    const double bodyPoleEclipticLongitude =   0.0*orsa::degToRad();
    
    // #warning OFF-North (km) and OFF-East (km) should be arguments!
    
#warning print obliquity...
    
#warning in most cases, the value of bodyRadius is not essential, only the orientation of the crater matters
    
    ORSA_DEBUG("thermal inertia: %g (SI)",thermalInertia());
    
#warning need to double-check these vectors when working at latitudes below equator
    
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
    // const orsa::Vector local_crater_center = 
    
    
    const size_t numSlices=400;
    History history;
    const size_t NS=100000;
    const size_t days=70;
    const double hdist = orsa::FromUnits(3.0,orsa::Unit::AU); // heliocentric distance
    ORSA_DEBUG("big-theta: %g",theta(hdist));
    std::vector<double> Fs;
    Fs.resize(NS);
    for (size_t p=0; p<NS; ++p) {
#warning restore this one!
        // Fs[p] = solar()/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2)*std::max(cos(days*orsa::twopi()*(double)p/(double)NS),0.0);
        
        // test
        /* Fs[p] = solar()/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2) *
           std::max(cos(days*orsa::twopi()*(double)p/(double)NS),0.0) *
           std::max(0.5+0.7*cos(orsa::twopi()*(double)p/(double)NS),0.0);
           if (p%99==0) Fs[p] = 0.0;
        */
        
        // TEST!
        double proj = std::max(cos(days*orsa::twopi()*(double)p/(double)NS),0.0);
        if (proj < 0.5) proj = 0.0;
        Fs[p] = solar()/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2)*proj;
        
    }
    
    /* for (size_t p=0; p<NS; ++p) {
       ORSA_DEBUG("Fs[%06i] = %g",p,Fs[p]);
       }
    */
    
    const double dx = 0.2*skinDepth(); // or a fraction of skinDepth = ls
    const double dt = days*rotationPeriod()/NS;
    const unsigned int history_skip = 1;
    
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
