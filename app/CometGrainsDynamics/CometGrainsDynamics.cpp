#include "CometGrainsDynamics.h"

int main (int argc, char **argv) {
    
    // input
    const double r_comet = orsa::FromUnits(1.0,orsa::Unit::AU);
    const double nucleus_ax = orsa::FromUnits(5.0,orsa::Unit::KM);
    const double nucleus_ay = orsa::FromUnits(4.0,orsa::Unit::KM);
    const double nucleus_az = orsa::FromUnits(3.0,orsa::Unit::KM);
    const size_t gravity_degree = 2;
    const double comet_density = orsa::FromUnits(orsa::FromUnits(1.0,orsa::Unit::GRAM),orsa::Unit::CM,-3);
    const double rotation_period = orsa::FromUnits(3.0,orsa::Unit::HOUR);
    const double pole_ecliptic_longitude =  0.0*orsa::degToRad();
    const double pole_ecliptic_latitude  = 90.0*orsa::degToRad();
    const double min_escape_velocity_factor = 0.5;
    const double max_escape_velocity_factor = 1.5;
    const double min_latitude = -90.0*orsa::degToRad();
    const double max_latitude = +90.0*orsa::degToRad();
    const double min_vertical_angle =  0.0*orsa::degToRad();
    const double max_vertical_angle = 45.0*orsa::degToRad();
    const double min_beta = 1.0e-6;
    const double max_beta = 1.0e-1;
    const int max_time_days = 100;
    
    const orsa::Time t0 = orsa::Time(0);
    const orsa::Time max_time(max_time_days,0,0,0,0);
    
    const double nucleus_volume = 4.0*orsa::pi()*nucleus_ax*nucleus_ay*nucleus_az/3.0;
    const double nucleus_mass = comet_density*nucleus_volume; 
    const double Hill_radius = orsa::HillRadius(r_comet,nucleus_mass,orsaSolarSystem::Data::MSun());
    const double omega = orsa::twopi()/rotation_period;
    
    orsa::Debug::instance()->initTimer();
    
    osg::ref_ptr<orsa::Body> sun = new orsa::Body;
    {
        sun->setName("sun");
        IBPS ibps;
        ibps.time = t0;
        ibps.inertial = new orsa::PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MSun());
        ibps.translational = new orsa::DynamicTranslationalBodyProperty;
        ibps.translational->setPosition(orsa::Vector(0,0,0));
        ibps.translational->setVelocity(orsa::Vector(0,0,0));
        sun->setInitialConditions(ibps);
    }
    
    osg::ref_ptr<orsa::Body> nucleus = new orsa::Body;
    osg::ref_ptr<orsa::Shape> nucleus_shape =
        new orsa::EllipsoidShape(nucleus_ax,nucleus_ay,nucleus_az);
    const orsa::Vector nucleus_r0 = orsa::Vector(r_comet,0,0);
    const orsa::Vector nucleus_v0 = orsa::Vector(0,sqrt(orsaSolarSystem::Data::GMSun()/r_comet),0); // circular orbit approximation, to keep the Hill sphere radius constant
    {
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
        ibps.rotational = new orsaSolarSystem::ConstantZRotationEcliptic_RotationalBodyProperty(t0,
                                                                                                0.0,
                                                                                                omega,
                                                                                                pole_ecliptic_longitude,
                                                                                                pole_ecliptic_latitude);
        nucleus->setInitialConditions(ibps);
    }

    size_t iter=0;
    while (1) {
        
        // loop on grains
        
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
        const orsa::Vector r0 = intersectionPoint;
        const orsa::Vector n0 = normal;
        
        const orsa::Vector u_rot =
            orsa::externalProduct(orsa::Vector(0,0,1),n0).normalized();
        const orsa::Vector u_pol =
            orsa::externalProduct(n0,u_rot).normalized();
        
        // horizontal direction of ejection of grain, measured from the direction of rotation
        const double phi = orsa::twopi()*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
        double s_phi, c_phi;
        sincos(phi,&s_phi,&c_phi);
        // vertical angle
        const double theta = min_vertical_angle + (max_vertical_angle-min_vertical_angle)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
        double s_theta, c_theta;
        sincos(theta,&s_theta,&c_theta);
        
        // not including rotation yet
#warning escape velocity approximate for points within the bounding sphere of the body
        const double escape_velocity = sqrt(2*orsa::Unit::G()*nucleus_mass/r0.length());
        const double ejection_velocity =
            escape_velocity*(min_escape_velocity_factor + (max_escape_velocity_factor-min_escape_velocity_factor)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform());
        
        // set velocity vector, including effect of nucleus rotation
        const orsa::Vector v0_rotational_component =
            orsa::externalProduct(orsa::Vector(0,0,omega),r0);
        const orsa::Vector v0 =
            n0*ejection_velocity*c_theta +
            u_rot*ejection_velocity*s_theta*c_phi +
            u_pol*ejection_velocity*s_theta*s_phi +
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
        
        
        osg::ref_ptr<orsa::Body> grain = new orsa::Body;
        {
            grain->setName("grain");
            IBPS ibps;
            ibps.time = t0;
            ibps.inertial = new orsa::PointLikeConstantInertialBodyProperty(0.0);
            ibps.translational = new orsa::DynamicTranslationalBodyProperty;
            ibps.translational->setPosition(r0+nucleus_r0);
            ibps.translational->setVelocity(v0+nucleus_v0);
            grain->beta = exp(log(min_beta) + (log(max_beta)-log(min_beta))*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform());
            grain->betaSun = sun.get();
            grain->setInitialConditions(ibps);
        }
        
        osg::ref_ptr<orsa::BodyGroup> bg = new BodyGroup;
        
        bg->addBody(sun);
        bg->addBody(nucleus);
        bg->addBody(grain);
        
        const double exo_radius = r_comet*sqrt(nucleus_mass/(grain->beta*orsaSolarSystem::Data::MSun()));
        const double bound_radius = std::min(Hill_radius,exo_radius);
        
        osg::ref_ptr<CGDIntegrator> integrator = new CGDIntegrator(grain.get(),nucleus.get(),bound_radius);
        integrator->integrate(bg.get(),
                              t0,
                              max_time,
                              orsa::Time(0,0,5,0,0));
        
        orsa::Time common_start_time, common_stop_time;
        const bool goodCommonInterval = bg->getCommonInterval(common_start_time,common_stop_time,false);
        if (!goodCommonInterval) {
            ORSA_DEBUG("problems...");
        } else {

            double final_distance;
            double lon_impact = -999*orsa::degToRad();
            double lat_impact = -999*orsa::degToRad();
            {
                const orsa::Time t = common_stop_time;
                orsa::Vector r,v;
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
            }

            // from beta to grain size
            const double    L = orsa::FromUnits(orsa::FromUnits(orsa::FromUnits(3.839e26,orsa::Unit::KG),orsa::Unit::METER,2),orsa::Unit::SECOND,-3); //
            const double    G = orsa::Unit::G();
            const double MSun = orsaSolarSystem::Data::MSun();
            const double    c = orsa::Unit::c();
            const double  Qpr = 1.0;
            //
            const double    rho_grain = comet_density; // or different?
            // Burns, Lamy, Soter 1979, Eq. (19)
            const double grain_radius = (3*L)/(16*orsa::pi()*G*MSun*c) * (Qpr)/(grain->beta*rho_grain);
            
            FILE * fp = fopen("CGD.out","a");
            gmp_fprintf(fp,"%g %g %g %g %.3e %.3e %.3e %.3e %g %g %g %g %7.3f %+7.3f %7.3f %7.3f %.3f %.3f %.3f %.3f %.3e %.3e %7.3f %.3e %.3e %i %8.3f %+8.3f\n",
                        orsa::FromUnits(r_comet,orsa::Unit::AU,-1),
                        orsa::FromUnits(nucleus_ax,orsa::Unit::KM,-1),
                        orsa::FromUnits(nucleus_ay,orsa::Unit::KM,-1),
                        orsa::FromUnits(nucleus_az,orsa::Unit::KM,-1),
                        orsa::FromUnits(nucleus_mass,orsa::Unit::KG,-1),
                        orsa::FromUnits(Hill_radius,orsa::Unit::KM,-1),
                        orsa::FromUnits(exo_radius,orsa::Unit::KM,-1),
                        orsa::FromUnits(bound_radius,orsa::Unit::KM,-1),
                        orsa::FromUnits(orsa::FromUnits(comet_density,orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
                        orsa::FromUnits(rotation_period,orsa::Unit::HOUR,-1),
                        pole_ecliptic_longitude*orsa::radToDeg(),
                        pole_ecliptic_latitude*orsa::radToDeg(),
                        lon*orsa::radToDeg(),
                        lat*orsa::radToDeg(),
                        phi*orsa::radToDeg(),
                        theta*orsa::radToDeg(),
                        orsa::FromUnits(orsa::FromUnits(escape_velocity,orsa::Unit::METER,-1),orsa::Unit::SECOND),
                        orsa::FromUnits(orsa::FromUnits(ejection_velocity,orsa::Unit::METER,-1),orsa::Unit::SECOND),
                        orsa::FromUnits(orsa::FromUnits(v0_rotational_component.length(),orsa::Unit::METER,-1),orsa::Unit::SECOND),
                        orsa::FromUnits(orsa::FromUnits(v0.length(),orsa::Unit::METER,-1),orsa::Unit::SECOND),
                        (*grain->beta),
                        orsa::FromUnits(grain_radius,orsa::Unit::CM,-1),
                        orsa::FromUnits(common_stop_time.get_d(),orsa::Unit::DAY,-1),
                        orsa::FromUnits((*integrator->max_distance),orsa::Unit::KM,-1),
                        orsa::FromUnits(final_distance,orsa::Unit::KM,-1),
                        integrator->outcome,
                        lon_impact*orsa::radToDeg(),
                        lat_impact*orsa::radToDeg());
            fclose (fp);
            
        }
        
        ++iter;
    }
    
    return 0;
}
