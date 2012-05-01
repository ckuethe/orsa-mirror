#include "CometGrainsDynamics.h"

int main (int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    ORSA_DEBUG("PID: %i",getpid());
    
    if (argc != 2) {
        ORSA_DEBUG("Usage: %s <runID>",argv[0]);
        exit(0);
    }
    
    const int runID = atoi(argv[1]);

    char filename[4096];
    sprintf(filename,"CGD_%02i.out",runID);
    const std::string filename_CGD = filename;
    sprintf(filename,"colden_%02i.out",runID);
    const std::string filename_colden = filename;
    
    // #warning comment out Random Seed in production
    
    // set randomSeed for testing purposes only
    // orsa::GlobalRNG::randomSeed = 1376174123;
    // orsa::GlobalRNG::randomSeed = 1119056643;
    // orsa::GlobalRNG::randomSeed = -128300218;
    // orsa::GlobalRNG::randomSeed = -124766705;
    // orsa::GlobalRNG::randomSeed = 1617326819;
    // orsa::GlobalRNG::randomSeed = 555;
    // orsa::GlobalRNG::randomSeed = 869523156;
    // orsa::GlobalRNG::randomSeed = -1316703212;
    // orsa::GlobalRNG::randomSeed = 1960578200;
    // orsa::GlobalRNG::randomSeed = 97180168;
    
    // NOTE: two alternative mechanisms for ejection velocity
    // 1) sampling distribution= rotational component + ejection velocity model (no gas drag)
    // 2) start with v=(rotational component only) and then gas drag increases it
    //
    // all depends on the gas_drag_coefficient value
    
#warning TODO: consider using fractal density for grains?
#warning TODO: mascons gravity
#warning TODO: fraction of active nucleus (list of lat-lon ranges) and correlate with total production rate
#warning TODO: fraction of grain as ice, fraction as dust
#warning TODO: 
    
    // input

    // 103P/Hartley 2
    const double comet_orbit_q = orsa::FromUnits(1.058690085281137,orsa::Unit::AU);
    const double comet_orbit_e = .6951452964967095;
    const double comet_orbit_i = 13.61716956119923*orsa::degToRad();
    const double comet_orbit_node = 219.7626609177958*orsa::degToRad();
    const double comet_orbit_peri = 181.1954811299036*orsa::degToRad();
    const orsa::Time comet_orbit_Tp = orsaSolarSystem::gregorTime(2010,10,28.25696720); // 2010-Oct-28.25696720
    const orsa::Time comet_orbit_epoch = orsaSolarSystem::gregorTime(2010,9,17.0); // comet_orbit_Tp; // orsaSolarSystem::gregorTime(2010,1,1);
    //
    const double nucleus_ax = orsa::FromUnits(1.00,orsa::Unit::KM);
    const double nucleus_ay = orsa::FromUnits(0.45,orsa::Unit::KM);
    const double nucleus_az = orsa::FromUnits(0.45,orsa::Unit::KM);
    const size_t gravity_degree = 2;
    const double comet_density = orsa::FromUnits(orsa::FromUnits(0.22,orsa::Unit::GRAM),orsa::Unit::CM,-3);
    const double grain_density = orsa::FromUnits(orsa::FromUnits(0.50,orsa::Unit::GRAM),orsa::Unit::CM,-3);
    const double rotation_period = orsa::FromUnits(18.34,orsa::Unit::HOUR);
    const double pole_phi_Tp = 0.0*orsa::degToRad(); // rotation angle at time Tp
    const double pole_ecliptic_longitude = +69.0*orsa::degToRad();
    const double pole_ecliptic_latitude  = +34.0*orsa::degToRad();
    const double min_latitude = -90.0*orsa::degToRad();
    const double max_latitude = +90.0*orsa::degToRad();
    const double min_grain_radius = orsa::FromUnits(0.000001,orsa::Unit::METER);
    const double max_grain_radius = orsa::FromUnits(1.000000,orsa::Unit::METER);    
    const int min_time_seconds =  60; // grains flying less than this time are not included
    const int max_time_days    = 200; // 100;
    const size_t pow_10_max_distance = 9;
    
    // gas (drag) coefficients
    const double gas_production_rate_at_1AU = orsa::FromUnits(1.0e28,orsa::Unit::SECOND,-1); // molecules/second
    const double gas_velocity_at_1AU = orsa::FromUnits(orsa::FromUnits(0.5,orsa::Unit::KM),orsa::Unit::SECOND,-1);
    const double gas_molar_mass = 18; // 18 for H20
    const double gas_drag_coefficient = 2.00; // Cd nominal: 0.40 (OR 2.00 ??)
    
    // molecules per unit area per unit time
#warning EYE ON THIS!!! (zero?)
#warning FACTOR for NON-pure ICE...
    const double grain_sublimation_rate = orsa::FromUnits(orsa::FromUnits(1.0e17,orsa::Unit::CM,-2),orsa::Unit::SECOND,-1);
    const double grain_sublimation_molecule_mass = orsa::FromUnits(gas_molar_mass*1.66e-27,orsa::Unit::KG); // conversion from molar
    
#warning drag coefficient Cd should be close to 2.0 when the grain size is close to the free mean path
    
    // const orsa::Time t_snapshot = comet_orbit_Tp - orsa::Time(60,0,0,0,0);
    const orsa::Time t_snapshot = orsaSolarSystem::gregorTime(2010,11,4.60); // comet_orbit_Tp; // + orsa::Time(30,0,0,0,0);
    
    const double nucleus_volume = 4.0*orsa::pi()*nucleus_ax*nucleus_ay*nucleus_az/3.0;
    const double nucleus_mass = comet_density*nucleus_volume; 
    const double omega = orsa::twopi()/rotation_period;
    
    orsaSPICE::SPICE::instance()->setDefaultObserver("SSB");
    orsaSPICE::SPICE::instance()->loadKernel("de405.bsp");
    
    // comet orbit
    orsaSolarSystem::OrbitWithEpoch comet_orbit;
    comet_orbit.mu = orsaSolarSystem::Data::GMSun();
    comet_orbit.e = comet_orbit_e;
    comet_orbit.a = comet_orbit_q/(1.0-comet_orbit.e);
    comet_orbit.i = comet_orbit_i;
    comet_orbit.omega_node       = comet_orbit_node;
    comet_orbit.omega_pericenter = comet_orbit_peri;
    comet_orbit.M                = twopi()*(comet_orbit_epoch-comet_orbit_Tp).get_d()/comet_orbit.period();
    comet_orbit.epoch = comet_orbit_epoch;
    
    osg::ref_ptr<orsa::BodyGroup> bg = new BodyGroup;
    
    osg::ref_ptr<orsa::Body> grain = new orsa::Body;
    IBPS grain_ibps;
    osg::ref_ptr<GrainUpdateIBPS> grainUpdateIBPS = new GrainUpdateIBPS;
    grain_ibps.updateIBPS = grainUpdateIBPS;
    
    size_t iter=0;
    while (iter < 1250000) {
        // while (1) {
        
        // start integration up to max_time_days before t_snapshot
        // const orsa::Time t0 = t_snapshot - orsa::Time(max_time_days,0,0,0,0)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
        // const orsa::Time t0 = t_snapshot - orsa::Time((max_time_days*86400)*(1000000*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform()));
        // #warning maybe use log scale for interval sampling!?
        // #warning have a minimumum here too? as it is, the minimum is 1 mu-sec...
        const orsa::Time t0 = t_snapshot - orsa::Time(exp(log(min_time_seconds*1e6) + (log(max_time_days*86400.0*1.0e6)-log(min_time_seconds*1e6))*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform()));
        
        // ORSA_DEBUG("t0: %20.12e",t0.get_d());
        
        // osg::ref_ptr<orsa::BodyGroup> bg = new orsa::BodyGroup;
        bg->clear();
        
        osg::ref_ptr<Body> sun = new orsa::Body;
        {
            sun->setName("SUN");
            sun->isLightSource = true;
            orsaSPICE::SpiceBodyTranslationalCallback * sbtc = new orsaSPICE::SpiceBodyTranslationalCallback(sun->getName());
            orsa::IBPS ibps;
            ibps.inertial = new PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MSun());
            ibps.translational = sbtc;
            sun->setInitialConditions(ibps);
        }
        bg->addBody(sun.get());
        
        osg::ref_ptr<Body> earth = new orsa::Body;
        {
            earth->setName("EARTH");
            orsaSPICE::SpiceBodyTranslationalCallback * sbtc = new orsaSPICE::SpiceBodyTranslationalCallback(earth->getName());
            orsa::IBPS ibps;
            ibps.inertial = new PointLikeConstantInertialBodyProperty(0.0); // no need for gravity perturbation from Earth, just need Earth position to project observations
            ibps.translational = sbtc;
            earth->setInitialConditions(ibps);
        }
        bg->addBody(earth);        
        
        osg::ref_ptr<orsa::Body> nucleus = new orsa::Body;
        osg::ref_ptr<orsa::EllipsoidShape> nucleus_shape;
        orsa::Vector nucleus_r0;
        orsa::Vector nucleus_v0;
        {
            nucleus_shape = new orsa::EllipsoidShape(nucleus_ax,nucleus_ay,nucleus_az);
            nucleus_shape->closestVertexEpsilonRelative = 1.0e-3;
            
            orsa::Vector rSun, vSun;
            if (!bg->getInterpolatedPosVel(rSun,vSun,sun.get(),t0)) {
                ORSA_DEBUG("problems...");
            }
            
            orsa::Vector rOrbit, vOrbit;
            {
                // linearly propagate the comet orbit to t0
                orsa::Orbit localOrbit = comet_orbit;
                localOrbit.M = fmod(comet_orbit.M + 
                                    orsa::twopi() * (t0-comet_orbit.epoch).get_d()/comet_orbit.period(),
                                    orsa::twopi());
                localOrbit.relativePosVel(rOrbit,vOrbit);
                
                rOrbit += rSun;
                vOrbit += vSun;
            }
            
            // 'export'
            nucleus_r0 = rOrbit;
            nucleus_v0 = vOrbit;
            
            nucleus->setName("nucleus");
            IBPS ibps;
            ibps.time = t0;
            // const double volume = 4.0*orsa::pi()*nucleus_ax*nucleus_ay*nucleus_az/3.0;
            // const double nucleus_mass = comet_density*volume;
            osg::ref_ptr<orsa::PaulMoment> pm = new orsa::PaulMoment(gravity_degree);
            orsa::EllipsoidExpansion(pm.get(),
                                     nucleus_ax,
                                     nucleus_ay,
                                     nucleus_az);
            if (0) {
                // test
                std::vector< std::vector<mpf_class> > C, S, norm_C, norm_S;
                std::vector<mpf_class> J;
                orsa::convert(C,
                              S,
                              norm_C,
                              norm_S,
                              J,
                              pm.get(),
                              orsa::FromUnits(10.0,orsa::Unit::KM),
                              true);
            }
            
            ibps.inertial = new orsa::ConstantInertialBodyProperty(nucleus_mass,
                                                                   nucleus_shape.get(),
                                                                   orsa::Vector(0,0,0),
                                                                   orsa::Matrix::identity(),
                                                                   orsa::Matrix::identity(),
                                                                   orsa::Matrix::identity(),
                                                                   pm.get());
            ibps.translational = new orsa::DynamicTranslationalBodyProperty;
            ibps.translational->setPosition(nucleus_r0);
            ibps.translational->setVelocity(nucleus_v0);
            // NOTE: initial angle will change again later...
            ibps.rotational = new orsaSolarSystem::ConstantZRotationEcliptic_RotationalBodyProperty(comet_orbit_Tp,
                                                                                                    pole_phi_Tp,
                                                                                                    omega,
                                                                                                    pole_ecliptic_longitude,
                                                                                                    pole_ecliptic_latitude);
            nucleus->setInitialConditions(ibps);
        }
        bg->addBody(nucleus);
        
        const double grain_initial_radius = exp(log(min_grain_radius) + (log(max_grain_radius)-log(min_grain_radius))*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform());
        const double grain_initial_beta   = GrainRadiusToBeta(grain_initial_radius,grain_density);
        
        // position of grain on the nucleus surface
        /* const double lon = orsa::twopi()*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
           WRONG: UNIFORM! const double lat = min_latitude  + (max_latitude-min_latitude)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
        */
        //
        double lon,lat;
        orsa::Vector r0,n0;
        size_t local_trials=0;
        const double min_cos_angle = cos(85.0*orsa::degToRad()); // less than 90 deg, not too close to 90 for efficiency
        if (min_cos_angle <= 0.0) {
            ORSA_DEBUG("problems...");
            exit(0);
        }
        const double max_r0 = std::max(nucleus_ax,std::max(nucleus_ay,nucleus_az));
        while (1) {
            ++local_trials;
            double x,y,z;
            orsa::GlobalRNG::instance()->rng()->gsl_ran_dir_3d(&x,&y,&z);
            const orsa::Vector u_ray(x,y,z);
            
            /* double s_lon, c_lon;
               sincos(lon,&s_lon,&c_lon);
               double s_lat, c_lat;
               sincos(lat,&s_lat,&c_lat);
               const orsa::Vector u_ray(c_lat*c_lon,
               c_lat*s_lon,
               s_lat);
            */
            orsa::Vector intersectionPoint;
            orsa::Vector normal;
            nucleus_shape->rayIntersection(intersectionPoint,
                                           normal,
                                           orsa::Vector(0,0,0),
                                           u_ray,
                                           false);
            // const orsa::Vector r0 = intersectionPoint;
            // just above the surface, to avoid roundoff problems
            r0 = intersectionPoint.normalized()*(intersectionPoint.length()+orsa::FromUnits(1.0,orsa::Unit::METER));
            n0 = normal;

            // this is needed to normalize for local projected surface area
            /* const double min_cos_angle = cos(85.0*orsa::degToRad()); // less than 90 deg, not too close to 90 for efficiency
               if (min_cos_angle <= 0.0) {
               ORSA_DEBUG("problems...");
               exit(0);
               }
            */
            // const double max_r0 = std::max(nucleus_ax,std::max(nucleus_ay,nucleus_az));
            const double cos_angle = std::max(min_cos_angle,r0.normalized()*n0.normalized());
            /* {
               static double max_angle = 0.0;
               max_angle = std::max(max_angle,acos(r0.normalized()*n0.normalized()));
               ORSA_DEBUG("lat: %g [deg]   lon: %g [deg] angle: %g [deg]   max: %g [deg]",
               lat*orsa::radToDeg(),
               lon*orsa::radToDeg(),
               acos(r0.normalized()*n0.normalized())*orsa::radToDeg(),
               max_angle*orsa::radToDeg());
               }
            */
            // const double rdm = (1.0/min_cos_angle)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
            const double rdm = (orsa::square(max_r0)/min_cos_angle)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
            if (rdm < r0.lengthSquared()/cos_angle) {
                // found (almost)
                lon = atan2(u_ray.getY(),u_ray.getX());
                lat = asin(u_ray.getZ());
                if ((lat >= min_latitude) && (lat <= max_latitude)) {
                    // ORSA_DEBUG("local_trials: %i",local_trials);
                    break; 
                }
            }
        }
        
        /* const orsa::Vector u_rot =
           orsa::externalProduct(orsa::Vector(0,0,1),n0).normalized();
           const orsa::Vector u_pol =
           orsa::externalProduct(n0,u_rot).normalized();
        */

        // not including rotation yet
#warning escape velocity approximate for points within the bounding sphere of the body
        const double escape_velocity = sqrt(2*orsa::Unit::G()*nucleus_mass/r0.length());
        
        // set velocity vector, including effect of nucleus rotation
        const orsa::Vector v0_rotational_component =
            orsa::externalProduct(orsa::Vector(0,0,omega),r0);
        const orsa::Vector v0 =
            v0_rotational_component;
        
        /* ORSA_DEBUG("lat: %g [deg]",lat*orsa::radToDeg());
           orsa::print(intersectionPoint);
           orsa::print(normal);
           orsa::print(u_rot);
           orsa::print(u_pol);
           ORSA_DEBUG("ejection velocity: %g [m/s] + rotational component: %g [m/s] = total velocity = %g [m/s]   [escape vel: %g [m/s]",
           ejection_velocity,v0_rotational_component.length(),v0.length(),escape_velocity);
           orsa::print(v0_rotational_component);
           orsa::print(v0);
        */
        
        // reset cV_l on ellipsoid shape
        // nucleus_shape->cV_l = 0.0;        
        
        // osg::ref_ptr<orsa::BodyGroup> bg = new BodyGroup;
        
        // osg::ref_ptr<orsa::Body> grain = new orsa::Body;
        // IBPS grain_ibps;
        {
            grain->setName("grain");
            // non-standard fields
            grainUpdateIBPS->reset();
            grainUpdateIBPS->bg = bg;
            grainUpdateIBPS->nucleus = nucleus;
            grainUpdateIBPS->grain = grain;
            grain_ibps.updateIBPS = grainUpdateIBPS;
            //
            grain_ibps.time = t0;
            grain_ibps.inertial = new GrainDynamicInertialBodyProperty(t0,
                                                                       grain_initial_radius,
                                                                       grain_density,
                                                                       grain_sublimation_rate,
                                                                       grain_sublimation_molecule_mass,
                                                                       grain.get());
            grain_ibps.translational = new orsa::DynamicTranslationalBodyProperty;
            const orsa::Matrix nucleus_l2g_t0 = orsa::localToGlobal(nucleus.get(),
                                                                    bg.get(),
                                                                    t0);
            grain_ibps.translational->setPosition(nucleus_r0+nucleus_l2g_t0*r0);
            grain_ibps.translational->setVelocity(nucleus_v0+nucleus_l2g_t0*v0);
            // grain->beta = grain_beta;
#warning need to keep updating beta...
            grain->beta = grain_initial_beta;
            grain->betaSun = sun.get();
            // gas drag
            // if (gas_drag_coefficient > 0.0) {
            {
                // this actuall also included sublimation force computation...
                grain->propulsion = new GasDrag(bg,
                                                sun,
                                                nucleus,
                                                grain,
                                                // grain->beta,
                                                grain_density,
                                                gas_production_rate_at_1AU,
                                                gas_velocity_at_1AU,
                                                gas_molar_mass,
                                                gas_drag_coefficient);
            }
            //
            grain->setInitialConditions(grain_ibps);
        }
        bg->addBody(grain);
        
        // gather some more initial conditions
        double initial_distance;
        double sun_initial_angle;
        double sun_initial_angle_360;
        double r_comet_t0;
        double Hill_radius;
        double exo_radius;
        double bound_radius;
        {
            const orsa::Time t = t0;
            orsa::Vector r,v;
            bg->getInterpolatedPosVel(r,
                                      v,
                                      sun,
                                      t);
            const orsa::Vector sun_r_global = r;
            const orsa::Vector sun_v_global = v;
            bg->getInterpolatedPosVel(r,
                                      v,
                                      nucleus,
                                      t);
            const orsa::Vector nucleus_r_global = r;
            const orsa::Vector nucleus_v_global = v;
            bg->getInterpolatedPosVel(r,
                                      v,
                                      grain,
                                      t);
            const orsa::Vector grain_r_global = r;
            const orsa::Vector grain_v_global = v;
            
            const orsa::Vector sun_r_relative_global = sun_r_global - nucleus_r_global;
            const orsa::Vector sun_v_relative_global = sun_v_global - nucleus_v_global;
            const orsa::Vector grain_r_relative_global = grain_r_global - nucleus_r_global;
            const orsa::Vector grain_v_relative_global = grain_v_global - nucleus_v_global;
            
            const orsa::Matrix g2l = orsa::globalToLocal(nucleus,bg,t);
            
            const orsa::Vector sun_r_relative_local = g2l*sun_r_relative_global;
            const orsa::Vector sun_v_relative_local = g2l*sun_v_relative_global;
            const orsa::Vector grain_r_relative_local = g2l*grain_r_relative_global;
            const orsa::Vector grain_v_relative_local = g2l*grain_v_relative_global;

            // local normal vectors
            const orsa::Vector u_sun = sun_r_relative_local.normalized();
            const orsa::Vector u_z   = orsa::Vector(0,0,1);
            const orsa::Vector u_ortho = orsa::externalProduct(u_z,u_sun).normalized();

            // orsa::print(u_sun);
            
            initial_distance = grain_r_relative_local.length();
            
            sun_initial_angle = acos((sun_r_global-nucleus_r_global).normalized()*grain_r_relative_global.normalized());
            
            sun_initial_angle_360 = atan2(grain_r_relative_local.normalized()*u_ortho,
                                          grain_r_relative_local.normalized()*u_sun);

            r_comet_t0 = sun_r_relative_global.length();
            
            Hill_radius  = orsa::HillRadius(r_comet_t0,nucleus_mass,orsaSolarSystem::Data::MSun());
            exo_radius   = r_comet_t0*sqrt(nucleus_mass/(grain->beta*orsaSolarSystem::Data::MSun()));
            bound_radius = std::min(Hill_radius,exo_radius);
            
            /* orsa::print((sun_r_global-nucleus_r_global).normalized()*grain_r_relative_global.normalized());
               orsa::print(sun_r_relative_local.normalized()*grain_r_relative_local.normalized());
               orsa::print(grain_r_relative_global.normalized());
               orsa::print(grain_r_relative_local.normalized());
               orsa::print(u_sun);
               orsa::print(u_ortho);
            */
        }
        
        osg::ref_ptr<CGDIntegrator> integrator = new CGDIntegrator(grain.get(),grain_initial_radius,grain_density,nucleus.get(),bound_radius,pow_10_max_distance,t0);
        // call singleStepDone once before starting, to perform initial checks
        orsa::Time dummy_time(0);
        integrator->singleStepDone(bg,t0,dummy_time,dummy_time);
        integrator->integrate(bg,
                              t0,
                              t_snapshot, // t0+orsa::Time(max_time_days,0,0,0,0),
                              orsa::Time(0,0,0,1,0));
        
        orsa::Time common_start_time, common_stop_time;
        const bool goodCommonInterval = bg->getCommonInterval(common_start_time,common_stop_time,false);
        if (!goodCommonInterval) {
            ORSA_DEBUG("problems...");
        } else {

            double final_distance;
            double lon_impact = -999*orsa::degToRad();
            double lat_impact = -999*orsa::degToRad();
            double sun_final_angle = -999*orsa::degToRad();
            {
                const orsa::Time t = common_stop_time;
                orsa::Vector r,v;
                bg->getInterpolatedPosVel(r,
                                          v,
                                          sun,
                                          t);
                const orsa::Vector sun_r_global = r;
                const orsa::Vector sun_v_global = v;
                bg->getInterpolatedPosVel(r,
                                          v,
                                          nucleus,
                                          t);
                const orsa::Vector nucleus_r_global = r;
                const orsa::Vector nucleus_v_global = v;
                bg->getInterpolatedPosVel(r,
                                          v,
                                          grain,
                                          t);
                const orsa::Vector grain_r_relative_global = r - nucleus_r_global;
                const orsa::Vector grain_v_relative_global = v - nucleus_v_global;
                const orsa::Matrix g2l = orsa::globalToLocal(nucleus,bg,t);
                const orsa::Vector grain_r_relative_local = g2l*grain_r_relative_global;
                const orsa::Vector grain_v_relative_local = g2l*grain_v_relative_global;
                
                // note: final_distance < body "radius" if impact happened (the impact is not resolved exactly at the edge)
                final_distance = grain_r_relative_local.length();
                
                if (integrator->outcome == CGDIntegrator::IMPACT) {
                    lon_impact = atan2(grain_r_relative_local.getY(),grain_r_relative_local.getX());
                    if (lon_impact < 0.0) lon_impact += orsa::twopi();
                    lat_impact = asin(grain_r_relative_local.getZ()/grain_r_relative_local.length());
                }
                
                sun_final_angle = acos((sun_r_global-nucleus_r_global).normalized()*grain_r_relative_global.normalized());
            }
            
            if (0) {

                FILE * fp = fopen(filename_CGD.c_str(),"a");
                gmp_fprintf(fp,"%.6f %g %g %g     %.3e %.3e %.3e %.3e %g     %g %g %g %g %g     %g %g %g %7.3f %+7.3f     %7.3f %7.3f %.3f %.3f %.3f     %.3f %.3e %.3e %10.6f %10.6f     %10.6f %10.6f %10.6f %10.6f %10.6f     %+.3e %+.3e %+.3e %+.3e %+.3e     %+.3e %.3e %.3e %.3e %i     %8.3f %+8.3f %+8.3f %+8.3f %+8.3f     %10.6f %i \n",
                            orsa::FromUnits(r_comet_t0,orsa::Unit::AU,-1),
                            orsa::FromUnits(nucleus_ax,orsa::Unit::KM,-1),
                            orsa::FromUnits(nucleus_ay,orsa::Unit::KM,-1),
                            orsa::FromUnits(nucleus_az,orsa::Unit::KM,-1),
                            /* 5 */ orsa::FromUnits(nucleus_mass,orsa::Unit::KG,-1),
                            orsa::FromUnits(Hill_radius,orsa::Unit::KM,-1),
                            orsa::FromUnits(exo_radius,orsa::Unit::KM,-1),
                            orsa::FromUnits(bound_radius,orsa::Unit::KM,-1),
                            orsa::FromUnits(orsa::FromUnits(comet_density,orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
                            /* 10 */ orsa::FromUnits(orsa::FromUnits(grain_density,orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
                            orsa::FromUnits(rotation_period,orsa::Unit::HOUR,-1),
                            pole_ecliptic_longitude*orsa::radToDeg(),
                            pole_ecliptic_latitude*orsa::radToDeg(),
                            gas_production_rate_at_1AU,
                            /* 15 */ gas_velocity_at_1AU,
                            gas_molar_mass,
                            gas_drag_coefficient,
                            lon*orsa::radToDeg(),
                            lat*orsa::radToDeg(),
                            /* 20 */ 0.0,
                            0.0,
                            orsa::FromUnits(orsa::FromUnits(escape_velocity,orsa::Unit::METER,-1),orsa::Unit::SECOND),
                            0.0, // orsa::FromUnits(orsa::FromUnits(ejection_velocity,orsa::Unit::METER,-1),orsa::Unit::SECOND),
                            orsa::FromUnits(orsa::FromUnits(v0_rotational_component.length(),orsa::Unit::METER,-1),orsa::Unit::SECOND),
                            /* 25 */ orsa::FromUnits(orsa::FromUnits(v0.length(),orsa::Unit::METER,-1),orsa::Unit::SECOND),
                            grain_initial_beta, // (*grain->beta),
                            orsa::FromUnits(grain_initial_radius,orsa::Unit::METER,-1),
                            orsa::FromUnits((common_stop_time-t0).get_d(),orsa::Unit::DAY,-1),
                            orsa::FromUnits((integrator->crossing_time[0]-t0).get_d(),orsa::Unit::DAY,-1),
                            /* 30 */ orsa::FromUnits((integrator->crossing_time[1]-t0).get_d(),orsa::Unit::DAY,-1),
                            orsa::FromUnits((integrator->crossing_time[2]-t0).get_d(),orsa::Unit::DAY,-1),
                            orsa::FromUnits((integrator->crossing_time[3]-t0).get_d(),orsa::Unit::DAY,-1),
                            0.0, // orsa::FromUnits((integrator->crossing_time[4]-t0).get_d(),orsa::Unit::DAY,-1),
                            0.0, // orsa::FromUnits((integrator->crossing_time[5]-t0).get_d(),orsa::Unit::DAY,-1),
                            /* 35 */ orsa::FromUnits(orsa::FromUnits(integrator->crossing_velocity[0],orsa::Unit::METER,-1),orsa::Unit::SECOND),
                            orsa::FromUnits(orsa::FromUnits(integrator->crossing_velocity[1],orsa::Unit::METER,-1),orsa::Unit::SECOND),
                            orsa::FromUnits(orsa::FromUnits(integrator->crossing_velocity[2],orsa::Unit::METER,-1),orsa::Unit::SECOND),
                            orsa::FromUnits(orsa::FromUnits(integrator->crossing_velocity[3],orsa::Unit::METER,-1),orsa::Unit::SECOND),
                            0.0, // orsa::FromUnits(orsa::FromUnits(integrator->crossing_velocity[4],orsa::Unit::METER,-1),orsa::Unit::SECOND),
                            /* 40 */ 0.0, // orsa::FromUnits(orsa::FromUnits(integrator->crossing_velocity[5],orsa::Unit::METER,-1),orsa::Unit::SECOND),
                            orsa::FromUnits(initial_distance,orsa::Unit::KM,-1),
                            orsa::FromUnits((*integrator->max_distance),orsa::Unit::KM,-1),
                            orsa::FromUnits(final_distance,orsa::Unit::KM,-1),
                            integrator->outcome,
                            /* 45 */ lon_impact*orsa::radToDeg(),
                            lat_impact*orsa::radToDeg(),
                            orsa::radToDeg()*sun_initial_angle,
                            orsa::radToDeg()*sun_initial_angle_360,
                            orsa::radToDeg()*sun_final_angle,
                            /* 50 */ orsa::FromUnits((t_snapshot-t0).get_d(),orsa::Unit::DAY,-1),
                            ((common_stop_time-t0) > (t_snapshot-t0)));
                fclose (fp);
            }
            
            // colden = column density file (including all the particles that do nor reach t_snapshot, to normalize production)
            {
                // init to 0 for grains that don't make it to t_snapshot
                double grainRadius = 0.0;
                double grainArea   = 0.0;
                //
                double grain_distance = 0.0;
                //
                double grain_R_X = 0.0;
                double grain_R_Y = 0.0;
                double grain_R_Z = 0.0;
                //
                double grain_R_sun = 0.0;
                double grain_R_orbit_pole  = 0.0;
                double grain_R_orbit_plane = 0.0;
                //
                double grain_R_orbit_velocity = 0.0;
                // grain_R_orbit_pole same as above;
                double grain_R_sunish = 0.0;
                //
                double grain_R_earth = 0.0;
                double grain_R_RA = 0.0;
                double grain_R_Dec = 0.0;
                //
                double grain_V = 0.0;
                //
                if ((common_stop_time-t0) >= (t_snapshot-t0)) {
                    
# warning much of this should be computed only once instead of in a loop...
                    
                    const orsa::Time t = t_snapshot;
                    orsa::Vector r,v;
                    bg->getInterpolatedPosVel(r,
                                              v,
                                              sun,
                                              t);
                    const orsa::Vector sun_r_global = r;
                    const orsa::Vector sun_v_global = v;
                    
                    bg->getInterpolatedPosVel(r,
                                              v,
                                              earth,
                                              t);
                    const orsa::Vector earth_r_global = r;
                    const orsa::Vector earth_v_global = v;
                    
                    bg->getInterpolatedPosVel(r,
                                              v,
                                              nucleus,
                                              t);
                    const orsa::Vector nucleus_r_global = r;
                    const orsa::Vector nucleus_v_global = v;
                    
                    // absolute ecliptic coordinate directions
                    const orsa::Vector u_X(1,0,0);
                    const orsa::Vector u_Y(0,1,0);
                    const orsa::Vector u_Z(0,0,1);

                    // sun, orbit pole, orbit plane
                    const orsa::Vector u_sun = (sun_r_global-nucleus_r_global).normalized();
                    const orsa::Vector u_tmp = (nucleus_v_global-sun_v_global).normalized();
                    const orsa::Vector u_orbit_pole  = orsa::externalProduct(u_tmp,u_sun).normalized();
                    const orsa::Vector u_orbit_plane = orsa::externalProduct(u_orbit_pole,u_sun).normalized();
                    // u_orbit_plane is in the general direction of comet velocity, but also orthogonal to orbit pole and sun direction

                    // similar, but velocity is always along velocity (not just orbit plane), but sun is not exact but "sunish"...
                    const orsa::Vector u_orbit_velocity = (nucleus_v_global-sun_v_global).normalized();
                    // u_orbit_pole is the same as above
                    const orsa::Vector u_sunish = orsa::externalProduct(u_orbit_pole,u_orbit_velocity).normalized();
                    
                    const orsa::Vector u_earth = (earth_r_global-nucleus_r_global).normalized();
                    const orsa::Vector u_tmp2  = orsaSolarSystem::equatorialToEcliptic()*orsa::Vector(0,0,1);
                    const orsa::Vector u_RA    = orsa::externalProduct(u_earth,u_tmp2).normalized();
                    const orsa::Vector u_Dec   = orsa::externalProduct(u_RA,u_earth).normalized();

                    static bool print=true;
                    if (print) {
                        
                        const orsa::Matrix nucleus_l2g_t = orsa::localToGlobal(nucleus.get(),
                                                                               bg.get(),
                                                                               t);
                        const orsa::Vector u_pole_global = nucleus_l2g_t*orsa::Vector(0,0,1);
                        ORSA_DEBUG("u_pole_global*u_X.............: %+.3f = %7.3f [deg]",u_pole_global*u_X,             acos(u_pole_global*u_X)*orsa::radToDeg());
                        ORSA_DEBUG("u_pole_global*u_Y.............: %+.3f = %7.3f [deg]",u_pole_global*u_Y,             acos(u_pole_global*u_Y)*orsa::radToDeg());
                        ORSA_DEBUG("u_pole_global*u_Z.............: %+.3f = %7.3f [deg]",u_pole_global*u_Z,             acos(u_pole_global*u_Z)*orsa::radToDeg());
                        ORSA_DEBUG("u_pole_global*u_sun...........: %+.3f = %7.3f [deg]",u_pole_global*u_sun,           acos(u_pole_global*u_sun)*orsa::radToDeg());
                        ORSA_DEBUG("u_pole_global*u_orbit_pole....: %+.3f = %7.3f [deg]",u_pole_global*u_orbit_pole,    acos(u_pole_global*u_orbit_pole)*orsa::radToDeg());
                        ORSA_DEBUG("u_pole_global*u_orbit_plane...: %+.3f = %7.3f [deg]",u_pole_global*u_orbit_plane,   acos(u_pole_global*u_orbit_plane)*orsa::radToDeg());
                        ORSA_DEBUG("u_pole_global*u_orbit_velocity: %+.3f = %7.3f [deg]",u_pole_global*u_orbit_velocity,acos(u_pole_global*u_orbit_velocity)*orsa::radToDeg());
                        ORSA_DEBUG("u_pole_global*u_sunish........: %+.3f = %7.3f [deg]",u_pole_global*u_sunish,        acos(u_pole_global*u_sunish)*orsa::radToDeg());
                        ORSA_DEBUG("u_pole_global*u_earth.........: %+.3f = %7.3f [deg]",u_pole_global*u_earth,         acos(u_pole_global*u_earth)*orsa::radToDeg());
                        ORSA_DEBUG("u_pole_global*u_RA............: %+.3f = %7.3f [deg]",u_pole_global*u_RA,            acos(u_pole_global*u_RA)*orsa::radToDeg());
                        ORSA_DEBUG("u_pole_global*u_Dec...........: %+.3f = %7.3f [deg]",u_pole_global*u_Dec,           acos(u_pole_global*u_Dec)*orsa::radToDeg());
                        print=false;
                    }
                    
                    bg->getInterpolatedPosVel(r,
                                              v,
                                              grain,
                                              t);
                    const orsa::Vector grain_r_relative_global = r - nucleus_r_global;
                    const orsa::Vector grain_v_relative_global = v - nucleus_v_global;
                    /* const orsa::Matrix g2l = orsa::globalToLocal(nucleus,bg,t);
                       const orsa::Vector grain_r_relative_local = g2l*grain_r_relative_global;
                       const orsa::Vector grain_v_relative_local = g2l*grain_v_relative_global;
                    */
                    const double obliquity = acos((orsa::localToGlobal(nucleus,bg,t)*orsa::Vector(0,0,1))*u_orbit_pole);
                    static bool print_obliquity=false;
                    if (!print_obliquity) {
                        ORSA_DEBUG("nucleus obliquity: %.3f [deg]",obliquity*orsa::radToDeg());
                        print_obliquity=true;
                    }
                    // orsa::print(orsa::localToGlobal(nucleus,bg,t));
                    
                    osg::ref_ptr <GrainDynamicInertialBodyProperty> inertial = dynamic_cast <GrainDynamicInertialBodyProperty*> (grain_ibps.inertial.get());
                    inertial->update(t_snapshot);
                    grainRadius = inertial->radius();
                    grainArea = orsa::pi()*orsa::square(grainRadius);
                    
                    grain_distance = grain_r_relative_global.length();
                    
                    grain_R_X = grain_r_relative_global*u_X;
                    grain_R_Y = grain_r_relative_global*u_Y;
                    grain_R_Z = grain_r_relative_global*u_Z;
                    
                    grain_R_sun         = grain_r_relative_global*u_sun;
                    grain_R_orbit_pole  = grain_r_relative_global*u_orbit_pole;
                    grain_R_orbit_plane = grain_r_relative_global*u_orbit_plane;

                    grain_R_orbit_velocity = grain_r_relative_global*u_orbit_velocity;
                    // grain_R_orbit_pole same as above
                    grain_R_sunish         = grain_r_relative_global*u_sunish;

                    grain_R_earth = grain_r_relative_global*u_earth;
                    grain_R_RA    = grain_r_relative_global*u_RA;
                    grain_R_Dec   = grain_r_relative_global*u_Dec;
                    
                    grain_V = grain_v_relative_global.length();
                    
                    // ORSA_DEBUG("rC: %i",inertial->referenceCount());
                }
                
                // ORSA_DEBUG("rC: %i",grain_ibps.inertial->referenceCount());
                
#warning note: some grains are behind the comet nucleus, so should not contribute to the column density...
                
#warning keep fields in sync with histo.cpp
                
                char line[4096];
                gmp_sprintf(line,"%.6f   %.5f %.5f   %.3e %.3e   %+8.3f %+7.3f   %.3e %.3e %.3e   %.3e   %+.3e %+.3e %+.3e   %+.3e %+.3e %+.3e   %+.3e %+.3e %+.3e   %+.3e %+.3e %+.3e   %.3e",
                            orsa::FromUnits(r_comet_t0,orsa::Unit::AU,-1),
                            //
                            orsaSolarSystem::timeToJulian(t0),
                            orsaSolarSystem::timeToJulian(t_snapshot),
                            //
                            orsa::FromUnits((t_snapshot-t0).get_d(),orsa::Unit::DAY,-1),
                            orsa::FromUnits((common_stop_time-t0).get_d(),orsa::Unit::DAY,-1),
                            //
                            lon*orsa::radToDeg(),
                            lat*orsa::radToDeg(),
                            //
                            orsa::FromUnits(grain_initial_radius,orsa::Unit::METER,-1),
                            orsa::FromUnits(grainRadius,orsa::Unit::METER,-1),
                            orsa::FromUnits(grainArea,orsa::Unit::METER,-2),
                            //
                            orsa::FromUnits(grain_distance,orsa::Unit::KM,-1),
                            //
                            orsa::FromUnits(grain_R_X,orsa::Unit::KM,-1),
                            orsa::FromUnits(grain_R_Y,orsa::Unit::KM,-1),
                            orsa::FromUnits(grain_R_Z,orsa::Unit::KM,-1),
                            //
                            orsa::FromUnits(grain_R_sun,orsa::Unit::KM,-1),
                            orsa::FromUnits(grain_R_orbit_pole,orsa::Unit::KM,-1),
                            orsa::FromUnits(grain_R_orbit_plane,orsa::Unit::KM,-1),
                            //
                            orsa::FromUnits(grain_R_orbit_velocity,orsa::Unit::KM,-1),
                            orsa::FromUnits(grain_R_orbit_pole,orsa::Unit::KM,-1),
                            orsa::FromUnits(grain_R_sunish,orsa::Unit::KM,-1),
                            //
                            orsa::FromUnits(grain_R_earth,orsa::Unit::KM,-1),
                            orsa::FromUnits(grain_R_RA,orsa::Unit::KM,-1),
                            orsa::FromUnits(grain_R_Dec,orsa::Unit::KM,-1),
                            //
                            orsa::FromUnits(orsa::FromUnits(grain_V,orsa::Unit::METER,-1),orsa::Unit::SECOND));
                
                // can add:
                // dr grain
                // dv grain
                // gas number density
                // gas mass density
                
                FILE * fp = fopen(filename_colden.c_str(),"a");
                gmp_fprintf(fp,"%s\n",line);
                fclose(fp);
            }
            
        }
        
        ++iter;
    }
    
    return 0;
}

bool GrainUpdateIBPS::update(const orsa::Time & t,
                             InertialBodyProperty * inertial,
                             TranslationalBodyProperty * translational,
                             RotationalBodyProperty * rotational) {
    
    // ORSA_DEBUG("called... t: %20.12e",t.get_d());
    
    // ORSA_DEBUG("set: %i",grain_r_relative_local_initial->var.isSet());
    
    orsa::UpdateIBPS::update(t,inertial,translational,rotational);
    
    // ORSA_DEBUG("bg: %x",bg.get());
    
    orsa::Vector r;
    bg->getInterpolatedPosition(r,
                                nucleus,
                                t);
    const orsa::Vector nucleus_r_global = r;
    
    // cannot use getInterpolatedPosition here...
    const orsa::Vector grain_r_global = translational->position();
    
    const orsa::Vector grain_r_relative_global = grain_r_global - nucleus_r_global;
    
    const orsa::Matrix g2l = orsa::globalToLocal(nucleus,bg,t);
    
    const orsa::Vector grain_r_relative_local = g2l*grain_r_relative_global;
    
    if (!grain_r_relative_local_initial->var.isSet()) {
        grain_r_relative_local_initial->var = grain_r_relative_local;
    } else {
#warning THIS SHOULD BE SMOOTH, and also DISTANCE THRESHOLD should be a PARAMETER!
        if ((grain_r_relative_local - grain_r_relative_local_initial->var).length() < orsa::FromUnits(100,orsa::Unit::METER)) {
            grain->beta = 0.0;
        } else {
            GrainDynamicInertialBodyProperty * inertial_cast =
                dynamic_cast<GrainDynamicInertialBodyProperty *> (inertial);
            grain->beta = GrainRadiusToBeta(inertial_cast->radius(),
                                            inertial_cast->_density);
        }
    }
    
    /* ORSA_DEBUG("called... t: %20.12e  distance: %g   beta: %g   now set: %i",
       t.get_d(),
       (grain_r_relative_local - grain_r_relative_local_initial->var).length(),
       (*grain->beta),
       grain_r_relative_local_initial->var.isSet());
    */
    
    return true;
}

