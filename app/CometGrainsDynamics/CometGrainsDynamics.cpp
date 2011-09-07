#include "CometGrainsDynamics.h"

int main (int argc, char **argv) {
    
    // input
    const double r_comet = orsa::FromUnits(1.0,orsa::Unit::AU);
    const double nucleus_ax = orsa::FromUnits(3.0,orsa::Unit::KM);
    const double nucleus_ay = orsa::FromUnits(2.5,orsa::Unit::KM);
    const double nucleus_az = orsa::FromUnits(2.0,orsa::Unit::KM);
    const double comet_density = orsa::FromUnits(orsa::FromUnits(1.0,orsa::Unit::GRAM),orsa::Unit::CM,-3);
    const double rotation_period = orsa::FromUnits(15.0,orsa::Unit::HOUR);
    const double pole_ecliptic_longitude =  0.0*orsa::degToRad();
    const double pole_ecliptic_latitude  = 90.0*orsa::degToRad();
    const double min_escape_velocity_factor = 0.8;
    const double max_escape_velocity_factor = 1.2;
    const double min_latitude = -90.0*orsa::degToRad();
    const double max_latitude = +90.0*orsa::degToRad();
    const double min_vertical_angle =  0.0*orsa::degToRad();
    const double max_vertical_angle = 60.0*orsa::degToRad();
    const double min_beta = 1.0e-6;
    const double max_beta = 1.0e-2;
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
        ibps.translational = new orsa::ConstantTranslationalBodyProperty(orsa::Vector(0,0,0));
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
        osg::ref_ptr<orsa::PaulMoment> pm = new orsa::PaulMoment(2);
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
    
    while (1) {
        
        // loop on grains
        
        osg::ref_ptr<orsa::Body> grain = new orsa::Body;
        {
            
            
            // position of grain on the nucleus surface
            const double lat = min_latitude  + (max_latitude-min_latitude)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
            const double lon = orsa::twopi()*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
            
            double s_lat, c_lat;
            sincos(lat,&s_lat,&c_lat);
            double s_lon, c_lon;
            sincos(lon,&s_lon,&c_lon);
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
            
            grain->setName("grain");
            IBPS ibps;
            ibps.time = t0;
            ibps.inertial = new orsa::PointLikeConstantInertialBodyProperty(0.0);
            ibps.translational = new orsa::DynamicTranslationalBodyProperty;
            ibps.translational->setPosition(r0+nucleus_r0);
            ibps.translational->setVelocity(v0+nucleus_v0);
            grain->beta = min_beta + (max_beta-min_beta)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();;
            grain->betaSun = sun.get();
            grain->setInitialConditions(ibps);
        }
        
        osg::ref_ptr<orsa::BodyGroup> bg = new BodyGroup;
        
        bg->addBody(sun);
        bg->addBody(nucleus);
        bg->addBody(grain);
        
        osg::ref_ptr<CGDIntegrator> integrator = new CGDIntegrator(grain.get(),nucleus.get());
        integrator->integrate(bg.get(),
                              t0,
                              max_time,
                              orsa::Time(0,0,5,0,0));
        
        orsa::Time common_start, common_stop;
        const bool goodCommonInterval = bg->getCommonInterval(common_start,common_stop,false);
        if (!goodCommonInterval) {
            ORSA_DEBUG("problems...");
        } else {
            print(common_stop);
            
        }
    }
    
    return 0;
}
