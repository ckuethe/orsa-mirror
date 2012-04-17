#include "CometGrainsDynamics.h"

int main (int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    // #warning comment out Random Seed in production
    
    // set randomSeed for testing purposes only
    // orsa::GlobalRNG::randomSeed = 1376174123;
    // orsa::GlobalRNG::randomSeed = 1119056643;
    // orsa::GlobalRNG::randomSeed = -128300218;
    // orsa::GlobalRNG::randomSeed = -124766705;
    // orsa::GlobalRNG::randomSeed = 1617326819;
    // orsa::GlobalRNG::randomSeed = 555;
    
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
    // const double r_comet = orsa::FromUnits(1.07,orsa::Unit::AU);
    const double comet_orbit_q = orsa::FromUnits(1.07,orsa::Unit::AU);
    const double comet_orbit_e = 0.2;
    const double comet_orbit_i = 25.0*orsa::degToRad();
    const double comet_orbit_node = 0.0;
    const double comet_orbit_peri = 0.0;
    const orsa::Time comet_orbit_Tp = orsaSolarSystem::gregorTime(2000,1,1.70);
    const orsa::Time comet_orbit_epoch = comet_orbit_Tp; // orsaSolarSystem::gregorTime(2010,1,1);
    //
    const double nucleus_ax = orsa::FromUnits(3.5,orsa::Unit::KM);
    const double nucleus_ay = orsa::FromUnits(2.4,orsa::Unit::KM);
    const double nucleus_az = orsa::FromUnits(2.2,orsa::Unit::KM);
    const size_t gravity_degree = 2;
    const double comet_density = orsa::FromUnits(orsa::FromUnits(0.4,orsa::Unit::GRAM),orsa::Unit::CM,-3);
    const double grain_density = orsa::FromUnits(orsa::FromUnits(0.5,orsa::Unit::GRAM),orsa::Unit::CM,-3);
    const double rotation_period = orsa::FromUnits(6.0,orsa::Unit::HOUR);
    const double pole_phi_Tp = 0.0*orsa::degToRad(); // rotation angle at time Tp
    const double pole_ecliptic_longitude =  0.0*orsa::degToRad();
    const double pole_ecliptic_latitude  = 90.0*orsa::degToRad();
    // const double min_ejection_velocity_constant = 0.5; // in the relation between beta and ejection velocity
    // const double max_ejection_velocity_constant = 1.5; // in the relation between beta and ejection velocity
    // const double ejection_velocity_beta_exponent = 0.5; // nominal: 0.5
    // const double ejection_velocity_radial_exponent = -0.5; // nominal: -0.5
    const double min_latitude = -90.0*orsa::degToRad();
    const double max_latitude = +90.0*orsa::degToRad();
    const double min_grain_radius = orsa::FromUnits(0.0001,orsa::Unit::METER);
    const double max_grain_radius = orsa::FromUnits(0.2000,orsa::Unit::METER);    
    const int max_time_days = 60; // 100;
    
    // gas drag coefficients
    const double gas_production_rate_at_1AU = orsa::FromUnits(1.0e28,orsa::Unit::SECOND,-1); // molecules/second
    const double gas_velocity_at_1AU = orsa::FromUnits(orsa::FromUnits(0.5,orsa::Unit::KM),orsa::Unit::SECOND,-1);
    const double gas_molar_mass = 18; // 18 for H20
    const double gas_drag_coefficient = 2.00; // Cd nominal: 0.40 (OR 2.00 ??)
    
    // molecules per unit area per unit time
#warning EYE ON THIS!!! (zero?)
    const double grain_sublimation_rate =  0.0; // orsa::FromUnits(orsa::FromUnits(1.0e17,orsa::Unit::CM,-2),orsa::Unit::SECOND,-1);
    const double grain_sublimation_molecule_mass = orsa::FromUnits(gas_molar_mass*1.66e-27,orsa::Unit::KG); // conversion from molar
    
#warning drag coefficient Cd should be close to 2.0 when the grain size is close to the free mean path
    
    // const orsa::Time t0 = orsa::Time(0);
    // const orsa::Time max_time(max_time_days,0,0,0,0);
    
    const orsa::Time t_snapshot = comet_orbit_Tp - orsa::Time(60,0,0,0,0);
    
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
    
    size_t iter=0;
    while (iter < 100000) {

        // start integration up to max_time_days before t_snapshot
        // const orsa::Time t0 = t_snapshot - orsa::Time(max_time_days,0,0,0,0)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
        const orsa::Time t0 = t_snapshot - orsa::Time((max_time_days*86400)*(1000000*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform()));
        
        // osg::ref_ptr<orsa::BodyGroup> bg = new BodyGroup;
        bg->clear();
        
        osg::ref_ptr<Body> sun = new Body;
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
        
        osg::ref_ptr<Body> earth = new Body;
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
            if (1) {
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
        const double lon = orsa::twopi()*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
        const double lat = min_latitude  + (max_latitude-min_latitude)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
        
        double s_lon, c_lon;
        sincos(lon,&s_lon,&c_lon);
        double s_lat, c_lat;
        sincos(lat,&s_lat,&c_lat);
        const orsa::Vector u_ray(c_lat*c_lon,
                                 c_lat*s_lon,
                                 s_lat);
        orsa::Vector intersectionPoint;
        orsa::Vector normal;
        nucleus_shape->rayIntersection(intersectionPoint,
                                       normal,
                                       orsa::Vector(0,0,0),
                                       u_ray,
                                       false);
        // const orsa::Vector r0 = intersectionPoint;
        // just above the surface, to avoid roundoff problems
        const orsa::Vector r0 = intersectionPoint*1.001;
        const orsa::Vector n0 = normal;
        
        const orsa::Vector u_rot =
            orsa::externalProduct(orsa::Vector(0,0,1),n0).normalized();
        const orsa::Vector u_pol =
            orsa::externalProduct(n0,u_rot).normalized();
        
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
        
        osg::ref_ptr<orsa::Body> grain = new orsa::Body;
        {
            grain->setName("grain");
            IBPS ibps;
            ibps.time = t0;
            ibps.inertial = new GrainDynamicInertialBodyProperty(t0,
                                                                 grain_initial_radius,
                                                                 grain_density,
                                                                 grain_sublimation_rate,
                                                                 grain_sublimation_molecule_mass,
                                                                 grain.get());
            ibps.translational = new orsa::DynamicTranslationalBodyProperty;
            const orsa::Matrix nucleus_l2g_t0 = orsa::localToGlobal(nucleus.get(),
                                                                    bg.get(),
                                                                    t0);
            ibps.translational->setPosition(nucleus_r0+nucleus_l2g_t0*r0);
            ibps.translational->setVelocity(nucleus_v0+nucleus_l2g_t0*v0);
            // grain->beta = grain_beta;
#warning need to keep updating beta...
            grain->beta = grain_initial_beta;
            grain->betaSun = sun.get();
            // gas drag
            if (gas_drag_coefficient > 0.0) {
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
            grain->setInitialConditions(ibps);
        }
        bg->addBody(grain);

        // NOT ANYMORE!
        /* 
           {
           // change initial angle of nucleus
           IBPS ibps = nucleus->getInitialConditions();        
           ibps.rotational = new orsaSolarSystem::ConstantZRotationEcliptic_RotationalBodyProperty(t0,
           orsa::twopi()*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform(),
           omega,
           pole_ecliptic_longitude,
           pole_ecliptic_latitude);
           nucleus->setInitialConditions(ibps);
           }
        */
        
        // compute at t_snapshot
        // const double Hill_radius = orsa::HillRadius( r_Comet   ,nucleus_mass,orsaSolarSystem::Data::MSun());
        // const double exo_radius = r_comet*sqrt(nucleus_mass/(grain->beta*orsaSolarSystem::Data::MSun()));
        // const double bound_radius = std::min(Hill_radius,exo_radius);
        
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
        
        osg::ref_ptr<CGDIntegrator> integrator = new CGDIntegrator(grain.get(),grain_initial_radius,grain_density,nucleus.get(),bound_radius,6,t0);
        // call singleStepDone once before starting, to perform initial checks
        orsa::Time dummy_time(0);
        integrator->singleStepDone(bg.get(),t0,dummy_time,dummy_time);
        integrator->integrate(bg.get(),
                              t0,
                              t0+orsa::Time(max_time_days,0,0,0,0),
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
            
            FILE * fp = fopen("CGD.out","a");
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
                        orsa::FromUnits((integrator->crossing_time[4]-t0).get_d(),orsa::Unit::DAY,-1),
                        orsa::FromUnits((integrator->crossing_time[5]-t0).get_d(),orsa::Unit::DAY,-1),
                        /* 35 */ orsa::FromUnits(orsa::FromUnits(integrator->crossing_velocity[0],orsa::Unit::METER,-1),orsa::Unit::SECOND),
                        orsa::FromUnits(orsa::FromUnits(integrator->crossing_velocity[1],orsa::Unit::METER,-1),orsa::Unit::SECOND),
                        orsa::FromUnits(orsa::FromUnits(integrator->crossing_velocity[2],orsa::Unit::METER,-1),orsa::Unit::SECOND),
                        orsa::FromUnits(orsa::FromUnits(integrator->crossing_velocity[3],orsa::Unit::METER,-1),orsa::Unit::SECOND),
                        orsa::FromUnits(orsa::FromUnits(integrator->crossing_velocity[4],orsa::Unit::METER,-1),orsa::Unit::SECOND),
                        /* 40 */ orsa::FromUnits(orsa::FromUnits(integrator->crossing_velocity[5],orsa::Unit::METER,-1),orsa::Unit::SECOND),
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
        
        ++iter;
    }
    
    return 0;
}

