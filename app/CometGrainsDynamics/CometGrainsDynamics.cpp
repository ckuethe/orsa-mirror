#include "CometGrainsDynamics.h"

int main (int argc, char **argv) {
    
    // input
    const double r_comet = orsa::FromUnits(1.0,orsa::Unit::AU);
    const double nucleus_ax = orsa::FromUnits(3.0,orsa::Unit::KM);
    const double nucleus_ay = orsa::FromUnits(2.5,orsa::Unit::KM);
    const double nucleus_az = orsa::FromUnits(2.0,orsa::Unit::KM);
    const double comet_density = orsa::FromUnits(orsa::FromUnits(1.0,orsa::Unit::GRAM),orsa::Unit::CM,-3);
    
    orsa::Debug::instance()->initTimer();
    
    osg::ref_ptr<Body> sun = new orsa::Body;
    {
        sun->setName("sun");
        IBPS ibps;
        ibps.time = orsa::Time(0);
        ibps.inertial = new orsa::PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MSun());
        ibps.translational = new orsa::ConstantTranslationalBodyProperty(orsa::Vector(0,0,0));
        sun->setInitialConditions(ibps);
    }
    
    osg::ref_ptr<Body> nucleus = new orsa::Body;
    {
        nucleus->setName("nucleus");
        IBPS ibps;
        ibps.time = orsa::Time(0);
        const double volume = 4.0*orsa::pi()*nucleus_ax*nucleus_ay*nucleus_az/3.0;
        const double nucleus_mass = comet_density*volume;
        osg::ref_ptr<orsa::PaulMoment> pm = new orsa::PaulMoment(2);
        orsa::EllipsoidExpansion(pm.get(),
                                 nucleus_ax,
                                 nucleus_ay,
                                 nucleus_az);
        ibps.inertial = new orsa::ConstantInertialBodyProperty(nucleus_mass,
                                                               new orsa::EllipsoidShape(nucleus_ax,nucleus_ay,nucleus_az),
                                                               orsa::Vector(0,0,0),
                                                               orsa::Matrix::identity(),
                                                               orsa::Matrix::identity(),
                                                               orsa::Matrix::identity(),
                                                               pm.get());
        ibps.translational = new orsa::ConstantTranslationalBodyProperty(orsa::Vector(r_comet,0,0));
        nucleus->setInitialConditions(ibps);
    }

    while (1) {
        // loop on grains
        
        osg::ref_ptr<orsa::BodyGroup> bg = new BodyGroup;
        
        bg->addBody(sun);
        bg->addBody(nucleus);
        // bg->addBody(grain);
        
        
    }
    
    return 0;
}
