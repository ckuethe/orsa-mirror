#ifndef COMET_GRAINS_DYNAMICS_H
#define COMET_GRAINS_DYNAMICS_H

#include <orsa/body.h>
#include <orsa/bodygroup.h>
#include <orsa/integrator_radau.h>
#include <orsa/integrator_leapfrog.h>
#include <orsa/orbit.h>
#include <orsa/paul.h>
#include <orsa/paulMoment.h>
#include <orsa/print.h>
#include <orsa/statistic.h>
#include <orsa/util.h>

#include <orsaSolarSystem/attitude.h>
#include <orsaSolarSystem/orbit.h>

#include <orsaSPICE/spice.h>
#include <orsaSPICE/spiceBodyTranslationalCallback.h>

// choose one depending on the shape file loaded
// #include "gaskell.h"
#include "gaskell_mod.h"

using namespace orsa;

// some global vars...

#warning use gas_plot_run = false in production...
static const bool gas_plot_run = true;

// Burns, Lamy, Soter 1979, Eq. (19)
double GrainBetaToRadius(const double & grainBeta,
                         const double & grainDensity,
                         const double & Qpr = 1.0) {
    static const double    L = orsa::FromUnits(orsa::FromUnits(orsa::FromUnits(3.839e26,orsa::Unit::KG),orsa::Unit::METER,2),orsa::Unit::SECOND,-3); //
    static const double    G = orsa::Unit::G();
    static const double MSun = orsaSolarSystem::Data::MSun();
    static const double    c = orsa::Unit::c();
    // static const double  Qpr = 1.0;
    return (3*L)/(16*orsa::pi()*G*MSun*c) * (Qpr)/(grainBeta*grainDensity);
}
double GrainRadiusToBeta(const double & grainRadius,
                         const double & grainDensity,
                         const double & Qpr = 1.0) {
    // actually, can call same function...
    return GrainBetaToRadius(grainRadius,grainDensity,Qpr);
}

/* double GrainRadius(const double & initialGrainRadius,
   const orsa::Time & t) {
   #warning COMPLETE THIS ONE!!!
   return initialGrainRadius;
   }
*/

// #error REPLACE  GrainIBPS  with  GrainUpdateIBPS 

// test, but this might be useful elsewhere... should improve and put in ORSA lib (if working like orsa::Cache<T> ...)
template <typename T> class RefCache : public osg::Referenced {
public:
    RefCache() : osg::Referenced() { }
protected:
    ~RefCache() { }
public:
    orsa::Cache<T> var;
};

class GrainUpdateIBPS : public orsa::UpdateIBPS {
public:
    GrainUpdateIBPS() : orsa::UpdateIBPS() {
        grain_r_relative_local_initial = new RefCache<orsa::Vector>;
    }
public:
    GrainUpdateIBPS(const GrainUpdateIBPS & grain_ibps) : orsa::UpdateIBPS() {
        bg = grain_ibps.bg;
        nucleus = grain_ibps.nucleus;
        grain = grain_ibps.grain;
        grain_r_relative_local_initial = grain_ibps.grain_r_relative_local_initial;
    }
protected:
    virtual ~GrainUpdateIBPS() { }
public:
    GrainUpdateIBPS & operator = (const GrainUpdateIBPS & grain_ibps) {
        bg = grain_ibps.bg;
        nucleus = grain_ibps.nucleus;
        grain = grain_ibps.grain;
        grain_r_relative_local_initial = grain_ibps.grain_r_relative_local_initial;
        return (*this);
    }
public:
    GrainUpdateIBPS * clone() const {
        return new GrainUpdateIBPS(*this);
    }
public:
    bool update(const orsa::Time & t,
                InertialBodyProperty * inertial,
                TranslationalBodyProperty * translational,
                RotationalBodyProperty * rotational);
public:
    void reset() const {
        grain_r_relative_local_initial->var.reset();
    }
public:
    osg::ref_ptr<orsa::BodyGroup> bg;
    osg::ref_ptr<orsa::Body> nucleus;
    osg::ref_ptr<orsa::Body> grain;
protected:
    // mutable orsa::Cache<orsa::Vector> grain_r_relative_local_initial;
    // osg::RefCache<orsa::Vector> grain_r_relative_local_initial;
    osg::ref_ptr< RefCache<orsa::Vector> > grain_r_relative_local_initial;
};

/* 
   class GrainIBPS : public orsa::IBPS {
   public:
   GrainIBPS() : orsa::IBPS() { }
   public:
   GrainIBPS(const GrainIBPS & grain_ibps) : orsa::IBPS(IBPS(grain_ibps)) {
   bg = grain_ibps.bg;
   nucleus = grain_ibps.nucleus;
   grain = grain_ibps.grain;
   grain_r_relative_local_initial = grain_ibps.grain_r_relative_local_initial;
   }
   public:
   virtual ~GrainIBPS() { }
   public:
   GrainIBPS & operator = (const GrainIBPS & grain_ibps) {
   (*this) = GrainIBPS(grain_ibps);
   bg = grain_ibps.bg;
   nucleus = grain_ibps.nucleus;
   grain = grain_ibps.grain;
   grain_r_relative_local_initial = grain_ibps.grain_r_relative_local_initial;
   }
   public:
   virtual bool update(const orsa::Time & t) {
   ORSA_DEBUG("called... t: %20.12e",t.get_d());
   orsa::IBPS::update(t);
   
   ORSA_DEBUG("bg: %x",bg.get());
   
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
   const orsa::Vector grain_r_global = r;
   const orsa::Vector grain_v_global = v;
   
   const orsa::Vector grain_r_relative_global = grain_r_global - nucleus_r_global;
   
   const orsa::Matrix g2l = orsa::globalToLocal(nucleus,bg,t);
   
   const orsa::Vector grain_r_relative_local = g2l*grain_r_relative_global;
   
   if (!grain_r_relative_local_initial.isSet()) {
   grain_r_relative_local_initial = grain_r_relative_local;
   } else {
   #warning THIS SHOULD BE SMOOTH, and also DISTANCE should be a PARAMETER!
   if ((grain_r_relative_local - grain_r_relative_local_initial).length() < orsa::FromUnits(100,orsa::Unit::METER)) {
   grain->beta = 0.0;
   }
   }
   
   return true;
   }
   public:
   osg::ref_ptr<orsa::BodyGroup> bg;
   osg::ref_ptr<orsa::Body> nucleus;
   osg::ref_ptr<orsa::Body> grain;
   protected:
   orsa::Cache<orsa::Vector> grain_r_relative_local_initial;
   };
*/

class GrainDynamicInertialBodyProperty : public InertialBodyProperty {
public:
    GrainDynamicInertialBodyProperty(const orsa::Time & t0,
                                     const double & initialRadius,
                                     const double & density,
                                     const double & sublimationRate_at_1AU, /* molecules per unit area per unit time */
                                     const double & moleculeMass, /* water molecule mass */
                                     const orsa::Body * grain,
                                     const orsa::Body * sun,
                                     const orsa::BodyGroup * bg) : 
        InertialBodyProperty(),
        _t0(t0),
        _initialRadius(initialRadius),
        _density(density),
        _dRadius_dt_at_1AU(sublimationRate_at_1AU*moleculeMass/(4.0*density)), /* check factors! i.e. the factor 4.0 due to energy balance, so it's sublimating from pi*Rg^2 instead of 4*pi*Rg^2 */
        _grain(grain),
        _sun(sun),
        _bg(bg) {
        // ORSA_DEBUG("dRdt: %g",_dRadius_dt);
        update_call_already_in_progress = false;
        _init();
    }
public:
    GrainDynamicInertialBodyProperty(const GrainDynamicInertialBodyProperty & ibp) : 
        InertialBodyProperty(ibp),
        _t0(ibp._t0),
        _initialRadius(ibp._initialRadius),
        _density(ibp._density),
        _dRadius_dt_at_1AU(ibp._dRadius_dt_at_1AU),
        _grain(ibp._grain),
        _sun(ibp._sun),
        _bg(ibp._bg) {
        update_call_already_in_progress = ibp.update_call_already_in_progress;
        _init();
    }
protected:
    virtual ~GrainDynamicInertialBodyProperty() { }
protected:
    void _init() {
        update(_t0);
    }
public:
    const orsa::Time _t0;
    const double _initialRadius;
    const double _density;
    const double _dRadius_dt_at_1AU; // radius shrinking rate
    const orsa::Body * _grain;
    const orsa::Body * _sun;
    const orsa::BodyGroup * _bg;
protected:
    double _radius; // updated by update(t) calls
    double _mass;   // updated by update(t) calls
public:
    double radius() const { return _radius; }
    double mass() const { return _mass; }
    const orsa::Shape * originalShape() const { return 0; }
    orsa::Vector centerOfMass() const { return orsa::Vector(0,0,0); }
    orsa::Matrix shapeToLocal() const { return orsa::Matrix::identity(); }
    orsa::Matrix localToShape() const { return orsa::Matrix::identity(); }
    orsa::Matrix inertiaMatrix() const { return orsa::Matrix::identity(); }
    const orsa::PaulMoment * paulMoment() const { return 0; }
    const orsa::MassCluster * massCluster() const { return 0; }
public:
    bool setMass(const double &) {
        ORSA_ERROR("this method should not have been called, please check your code.");
        // orsa::crash();
        return false;
    }
    bool setOriginalShape(const orsa::Shape *) {
        ORSA_ERROR("this method should not have been called, please check your code.");
        return false;
    }
    bool setCenterOfMass(const orsa::Vector &) {
        ORSA_ERROR("this method should not have been called, please check your code.");
        return false;
    }
    bool setShapeToLocal(const orsa::Matrix &) {
        ORSA_ERROR("this method should not have been called, please check your code.");
        return false;
    }
    bool setLocalToShape(const orsa::Matrix &) {
        ORSA_ERROR("this method should not have been called, please check your code.");
        return false;
    }
    bool setInertiaMatrix(const orsa::Matrix &) {
        ORSA_ERROR("this method should not have been called, please check your code.");
        return false;
    }
    bool setPaulMoment(const orsa::PaulMoment *) {
        ORSA_ERROR("this method should not have been called, please check your code.");
        return false;
    }
    bool setMassCluster(const orsa::MassCluster *) {
        ORSA_ERROR("this method should not have been called, please check your code.");
        return false;
    }
public:
    InertialBodyProperty * clone() const {
        return new GrainDynamicInertialBodyProperty(*this);
    }
public:
    BodyPropertyType type() const { return BP_DYNAMIC; }

protected:
    // hack to prevent recursive calls...
    static bool update_call_already_in_progress;
public:
    bool update(const orsa::Time & t) {
        
        /* ORSA_DEBUG("update_call_already_in_progress: %i   this: %x",
           update_call_already_in_progress,this);
        */
        
        if (update_call_already_in_progress) return true;
        update_call_already_in_progress = true;
        
#warning minimum radius here should be a parameter...
        const double min_radius = orsa::FromUnits(1.0e-6,orsa::Unit::METER);
        
        // orsa::print(t);
        
        const double dt = (t-_t0).get_d();
        
        // ORSA_DEBUG("dt: %g",dt);
        
        if (t==_t0) {
            _radius = std::max(min_radius,_initialRadius);
        } else {
            orsa::Vector rSun, vSun;
            if (!_bg->getInterpolatedPosVel(rSun,vSun,_sun,t)) {
                ORSA_DEBUG("problems...");
            }	
            orsa::Vector rGrain,vGrain;
            if (!_bg->getInterpolatedPosVel(rGrain,vGrain,_grain,t)) {
                ORSA_DEBUG("problems...");
            }
            /* 
               orsa::IBPS ibps;
               _bg->getIBPS(ibps,_grain,t);
               ibps.lock();
               const orsa::Vector rGrain = ibps.translational->position();
               ibps.unlock();
            */
            const orsa::Vector R_h = (rGrain-rSun);
            const double r_h = R_h.length();
            const double r_h_AU = orsa::FromUnits(r_h,orsa::Unit::AU,-1);
            
            _radius = std::max(min_radius,_initialRadius-_dRadius_dt_at_1AU*pow(r_h_AU,-2)*dt);
        }
        
        // ORSA_DEBUG("dt: %g  _radius: %g",dt,_radius);
        
#warning check this and IMPROVE... (slow sublimation at 10 microns)
#warning this is not correct for icy grains that have an initial radius BELOW 10 microns
#warning written as is, is wrong, especially for dust-only grains (not sublimating)
        /* if (_radius < orsa::FromUnits(1.0e-5,orsa::Unit::METER)) {
           _radius = orsa::FromUnits(1.0e-5,orsa::Unit::METER);
           }
        */
        
        _mass = 4.0/3.0*orsa::pi()*orsa::cube(_radius)*_density;
        
        // _grain->beta is now updated in GrainUpdateIBPS
        
        update_call_already_in_progress = false;
        
        return true;
    }
public:
    void lock() { }
    void unlock() { }
};

class GasDrag : public orsa::Propulsion {
public:
    GasDrag(orsa::BodyGroup  * bg_in,
            const orsa::Body * sun_in,
            const orsa::Body * comet_in,
            const orsa::Body * grain_in,
            // const double & grain_beta_in,
            const double & grain_density_in,
            const double & nucleus_water_production_rate_factor, // gas_production_rate_at_1AU,
            const double & gas_velocity_at_1AU,
            const double & gas_molar_mass, // i.e. 18 for H20
            const double & gas_drag_coefficient) :
        orsa::Propulsion(),
        bg(bg_in),
        sun(sun_in),
        comet(comet_in),
        grain(grain_in),
        // grainBeta(grain_beta_in),
        grainDensity(grain_density_in),
        // grainRadius(GrainBetaToRadius(grainBeta,grainDensity)),
        // grainArea(orsa::pi()*orsa::square(grainRadius)),
        eta_water(nucleus_water_production_rate_factor),
        Vgas_1AU(gas_velocity_at_1AU),
        Mgas(orsa::FromUnits(gas_molar_mass*1.66e-27,orsa::Unit::KG)), // conversion from molar
        Cd(gas_drag_coefficient),
        newton(orsa::FromUnits(orsa::FromUnits(orsa::FromUnits(1,orsa::Unit::KG),orsa::Unit::METER),orsa::Unit::SECOND,-2)) { }
protected:
    virtual ~GasDrag() { }
public:	
    orsa::Vector getThrust(const orsa::Time & t) const {
        
        static const double water_sublimation_rate_at_1AU = orsa::FromUnits(orsa::FromUnits(1.0e17,orsa::Unit::CM,-2),orsa::Unit::SECOND,-1);
        
        // ORSA_DEBUG("referenceCount(): %i",referenceCount());
        
#warning this would be good if force from sublimation was not computed too in this method
        // if (Cd == 0.0) return orsa::Vector(0,0,0);
        
        orsa::Vector rSun, vSun;
        if (!bg->getInterpolatedPosVel(rSun,vSun,sun,t)) {
            ORSA_DEBUG("problems...");
        }	
        
        orsa::Vector rComet, vComet;
        if (!bg->getInterpolatedPosVel(rComet,vComet,comet,t)) {
            ORSA_DEBUG("problems...");
        }
        
        orsa::Vector rGrain,vGrain;
        if (!bg->getInterpolatedPosVel(rGrain,vGrain,grain,t)) {
            ORSA_DEBUG("problems...");
        }	

        // triplet of orthogonal unit vectors: uS (towards Sun), uV (along comet velocity), uN (normal to comet orbit plane)
        const orsa::Vector uS = (rSun-rComet).normalized();
        const orsa::Vector u_tmp_z = orsa::externalProduct(uS,(vSun-vComet)).normalized();
        const orsa::Vector uV = orsa::externalProduct(u_tmp_z,uS).normalized();
        const orsa::Vector uN = orsa::externalProduct(uS,uV).normalized();
        
        /* ORSA_DEBUG("---- %g %g %g",uS*uV,uS*uN,uV*uN);
           orsa::print(uS);
           orsa::print(uV);
           orsa::print(uN);
        */
        
        const orsa::Matrix g2l = orsa::globalToLocal(comet,bg,t);
        const orsa::Matrix l2g = orsa::localToGlobal(comet,bg,t);
                
        orsa::IBPS ibps;
        if (!bg->getIBPS(ibps,comet,t)) {
            ORSA_DEBUG("problems at t:");
            orsa::print(t);
            orsa::crash();
        }
        /* 
           const orsa::EllipsoidShape * nucleus_shape =
           dynamic_cast<const orsa::EllipsoidShape *> (ibps.inertial->originalShape());
           double na, nb, nc;
           nucleus_shape->getABC(na,nb,nc);
           const double nucleus_max_radius = std::max(na,std::max(nb,nc));
        */
        const GaskellPlateModel * nucleus_shape =
            dynamic_cast<const GaskellPlateModel *> (ibps.inertial->originalShape());
        if (!nucleus_shape) {
            ORSA_DEBUG("problems...");
            exit(0);
        }
        const double nucleus_volume_equivalent_average_radius = cbrt(3*nucleus_shape->volume()/(4*orsa::pi()));
        /* ORSA_DEBUG("nucleus_volume_equivalent_average_radius = %g [km]   nucleus volume: %g [km^3]",
           orsa::FromUnits(nucleus_volume_equivalent_average_radius,orsa::Unit::KM,-1),
           orsa::FromUnits(nucleus_shape->volume(),orsa::Unit::KM,-3));
        */
        
        const orsa::TriShape::VertexVector & vv = nucleus_shape->getVertexVector();
        
#warning this would be good if force from sublimation was not computed too in this method
        /* 
           if (nucleus_shape->isInside(g2l*(rGrain-rComet))) {
           // ORSA_DEBUG("zero thrust for grain inside comet nucleus...");
           // return orsa::Vector(0,0,0);
           }
        */
        
        // switch to radial when difference is approximately smaller than 1 deg
        // orsa::Vector u_gas;
        /* #warning shoud replace this test with one more physically sound, comparing rotation period with gas expansion time to reach the distance
           if (r_c/nucleus_max_radius > 1000) {
           // radial
           u_gas = R_c.normalized();
           } else {
        */
        // orsa::Vector n0;
        /* 
           {
           // relative to comet
           // const orsa::Vector V_Gas_c   = v_gas_h * (rGrain-rComet).normalized();
           // modify V_gas_c to smoothly decrease near nucleus
           const orsa::Vector dr_g = rGrain-rComet;
           const orsa::Vector dr_l = g2l*dr_g;
           const orsa::Vector normal_l = nucleus_shape->_getVertexNormal(nucleus_shape->closestVertexIndex(dr_l));
           n0 = normal_l.normalized();
           const orsa::Vector normal_g = l2g*normal_l;
           u_gas = normal_g;
           }
        */
        
        const orsa::Vector R_h = (rComet-rSun);
        const double r_h = R_h.length();
        
        const double r_h_AU = orsa::FromUnits(r_h,orsa::Unit::AU,-1);
        
        const orsa::Vector R_c = (rGrain-rComet);
        const double r_c = R_c.length();
        
        const double water_sublimation_rate = water_sublimation_rate_at_1AU * pow(r_h_AU,-2);
        
        // const double dist_ratio = r_c / nucleus_shape->boundingRadius();
        const double dist_ratio = r_c / nucleus_volume_equivalent_average_radius;
        // gas velocity at r_h, relative to comet
        const double v_gas_h = Vgas_1AU * pow(r_h_AU,-0.5);
        const double v_gas_factor = dist_ratio/(1.0+dist_ratio); // goes from 0.5 near nucleus to 1.0 asymptotically
        const double v_gas_at_grain = v_gas_h*v_gas_factor;
        
        // ORSA_DEBUG("v_gas_at_grain: %g [km/s]",orsa::FromUnits(orsa::FromUnits(v_gas_at_grain,orsa::Unit::KM,-1),orsa::Unit::SECOND));
        
        // compute at class init?
        const double production_rate_per_unit_surface = eta_water*water_sublimation_rate;
        
#warning compute plate center and size once for all in Shape class?
        
        double effective_production_rate=0.0;
        orsa::Vector tmp_u_gas(0,0,0);
        double tmp_number_density_gas=0.0;
        // double max_test=0.0, mtt1=0.0, mtt2=0.0, mtt3=0.0, mtt4=0.0;
        for (size_t plt=0; plt<nucleus_shape->getFaceVector().size(); ++plt) {
            
            const orsa::TriShape::TriIndex plt_TriIndex = nucleus_shape->getFaceVector()[plt];
            
            const orsa::Vector & v_i = vv[plt_TriIndex.i()];
            const orsa::Vector & v_j = vv[plt_TriIndex.j()];
            const orsa::Vector & v_k = vv[plt_TriIndex.k()];
            
            const orsa::Vector & plt_center = (v_i+v_j+v_k)/3.0;
            
            const orsa::Vector & plt_normal = nucleus_shape->_getFaceNormal(plt);
            
            const orsa::Vector delta_l = (g2l*(rGrain-rComet)-plt_center);
            
            const double theta_sun = acos(g2l*(rSun-rComet).normalized() * plt_normal);
            
            // const double theta_sun_factor = std::max(0.0,cos(0.5*theta_sun));
            // const double theta_sun_factor = std::max(0.0,pow(cos(theta_sun),0.25)); // temperature for a low thermal inertia body goes as cos(theta_sun)^(1/4)
            const double theta_sun_factor = std::max(0.0,cos(theta_sun));
            
#warning compare factors all around in this loop and make sure they match
            
            effective_production_rate +=
                production_rate_per_unit_surface *
                nucleus_shape->_getFaceArea(plt) *
                theta_sun_factor;
            
            if (plt_normal*delta_l <= 0.0) continue;
            
            // ORSA_DEBUG("one good...");
            
            // RMS of each side lenght
            /* const double plt_RMS_scale = sqrt(((v_j-v_i).lengthSquared()+
               (v_k-v_i).lengthSquared()+
               (v_k-v_j).lengthSquared())/3.0);
            */

            // const double mx_factor = std::max(0.0,plt_normal*delta_l) / pow(delta_l.lengthSquared()+orsa::square(plt_RMS_scale),1.5);
            // const double mx_factor = std::max(0.0,plt_normal*delta_l) / pow(delta_l.lengthSquared()+nucleus_shape->_getFaceArea(plt),1.5);
            // const double mx_factor = std::max(0.0,plt_normal*delta_l) * nucleus_shape->_getFaceArea(plt) / pow(delta_l.lengthSquared()+nucleus_shape->_getFaceArea(plt),1.5);
            // const double mx_factor = 1.0 / (nucleus_shape->_getFaceArea(plt) + orsa::square(plt_normal*delta_l));
            // const double lp = (delta_l.length()-delta_l*plt_normal);
            const double lp = (delta_l-(delta_l*plt_normal)*plt_normal).length();
            const double lp_ratio  = lp/sqrt(nucleus_shape->_getFaceArea(plt)/orsa::pi());
            const double lp_factor = 1.0 / (1.0 + orsa::square(lp_ratio));
            const double lp_integral_normalization = 1.0/orsa::square(orsa::pi());
            /* const double mx_factor =
               lp_integral_normalization * lp_factor / (nucleus_shape->_getFaceArea(plt) + orsa::square(plt_normal*delta_l)); // ~ 1/r^2
            */
            const double mx_factor =
                lp_integral_normalization * lp_factor / (nucleus_shape->_getFaceArea(plt) + delta_l.lengthSquared()); // ~ 1/r^2
            
            
            
            
            
#warning fix cos theta sun
            
            const double plt_term = 
                production_rate_per_unit_surface *
                nucleus_shape->_getFaceArea(plt) *
                theta_sun_factor *
                mx_factor;
            
            // std::max(0.0,plt_normal*delta_l) / pow(orsa::square(plt_normal*delta_l)+orsa::square(plt_RMS_scale),1.5);
            // std::max(0.0,plt_normal*delta_l) / pow(orsa::square(plt_normal*delta_l)+1.0,1.5);
            // std::max(0.0,plt_normal*delta_l) / pow(delta_l.lengthSquared()+orsa::square(plt_RMS_scale),1.5);
            tmp_number_density_gas += plt_term;
            // tmp_u_gas += plt_term*delta_l.normalized();
            // tmp_u_gas += plt_term*plt_normal;
            // tmp_u_gas += plt_term*(lp_factor*plt_normal+(1.0-lp_factor)*delta_l.normalized());
            tmp_u_gas += plt_term*delta_l.normalized();
            
            /* if (plt_term > max_test) {
               max_test = plt_term;
               mtt1 = delta_l.length();
               mtt2 = lp_ratio;
               mtt3 = lp_factor;
               mtt4 = nucleus_shape->_getFaceArea(plt) * mx_factor;
               }   
            */
            
        }
        const double number_density_gas_at_grain = 1.0/(orsa::twopi()*v_gas_at_grain)*tmp_number_density_gas;
        // const orsa::Vector u_gas_l = 1.0/(orsa::twopi()*v_gas_at_grain*number_density_gas_at_grain)*tmp_u_gas.normalized();
        // ORSA_DEBUG("------------");
        const orsa::Vector u_gas_l = tmp_u_gas.normalized();
        const orsa::Vector v_gas_g = l2g*u_gas_l*v_gas_at_grain;
        // orsa::print(u_gas_l);
        // orsa::print(v_gas_g);
        
        /* ORSA_DEBUG("max_test = %10.5e mtt: %10.5f %10.5f %10.5f %10.5f",
           max_test,mtt1,mtt2,mtt3,mtt4);
        */
        
        const double rho = number_density_gas_at_grain * Mgas;
        
        // ORSA_DEBUG("effecitive production rate: %g [molecules/second]",orsa::FromUnits(effective_production_rate,orsa::Unit::SECOND));
        // ORSA_DEBUG("number density: %g",number_density_gas_at_grain);
        
        // ORSA_DEBUG("scp: %g",v_gas_g.normalized()*(rGrain-rComet).normalized());
        
        // gas velocity at r_h, relative to comet
        // const double v_gas_h = Vgas_1AU * pow(r_h_AU,-0.5);
        
        // production rate at r_h
        // const double Q_h = Q_1AU * pow(r_h_AU,-2);
        
        // number density at r_c for production rate Q_h
        //const double n = Q_h / (4*orsa::pi()*orsa::square(r_c)*v_gas_h);
        // optional: can multiply x cos(theta_sun) to account for gas only from lit side of comet
        // const double theta_sun = acos((rSun-rComet).normalized() * R_c.normalized());
        //         const double theta_sun = acos((rSun-rComet).normalized() * n0);
        // #warning MAKE SURE n0 is LOCAL or GLOBAL...
        
        
        // #warning restore this one!
        // const double theta_sun_factor = cos(0.5*theta_sun);
        // const double theta_sun_factor = 1.0;
        
        // const double theta_sun_factor = std::max(0.0,pow(cos(theta_sun),0.25)); // temperature for a low thermal inertia body goes as cos(theta_sun)^(1/4)
        // const double n = Q_h * theta_sun_factor / (4*orsa::pi()*orsa::square(r_c)*v_gas_h);
        
        // ORSA_DEBUG("theta_sun_factor: %g",theta_sun_factor);
        
        // rho = mass density = number density x molecular mass
        // const double rho = n * Mgas;
        
        orsa::IBPS grain_ibps;
        // GrainIBPS grain_ibps;
        bg->getIBPS(grain_ibps,grain,t);
        osg::ref_ptr <GrainDynamicInertialBodyProperty> inertial = dynamic_cast <GrainDynamicInertialBodyProperty*> (grain_ibps.inertial.get());
        const double grainRadius = inertial->radius();
        // if (grainRadius <= 0.0) return orsa::Vector(0,0,0);
        const double grainArea = orsa::pi()*orsa::square(grainRadius);
        
        // NOTE: gas direction is proportional to normal_g, which gets close to radial at large distances
        
        // const double dist_ratio = r_c / nucleus_max_radius;
        // const double dist_ratio = r_c / nucleus_volume_equivalent_average_radius;
        // const double v_Gas_factor = dist_ratio/(1.0+dist_ratio); // goes from 0.5 near nucleus to 1.0 asymptotically
        // const orsa::Vector V_Gas_c = v_Gas_factor * v_gas_h * u_gas;
        
        const orsa::Vector V_Grain_c = vGrain-vComet;
        // const double dV = (V_Grain_c - V_Gas_c)*(V_Gas_c.normalized());
        const double dV = (V_Grain_c - v_gas_g) * v_gas_g.normalized();
        const double sign = (dV>0) ? -1 : +1;
        const orsa::Vector thrust =
            sign*0.5*rho*dV*dV*Cd*grainArea*v_gas_g.normalized();
        
        /* ORSA_DEBUG("dV: %g m/s",orsa::FromUnits(orsa::FromUnits(dV,orsa::Unit::METER,-1),orsa::Unit::SECOND));
           orsa::print(V_Grain_c);
           orsa::print(v_gas_g);
        */
        
        // force due to sublimation? acc = mol*Vgas*Z / (Rg*rho_grain)
        orsa::Vector sublimationForce(0,0,0);
        if (0) {
#warning MUST pass the parameters as arguments...
#warning what is the multiplicative factor in front?
            const double sublimationForceFactor = 0.00;
            const double grain_sublimation_rate_at_1AU = orsa::FromUnits(orsa::FromUnits(1.0e17,orsa::Unit::CM,-2),orsa::Unit::SECOND,-1);
            sublimationForce =
                -uS*sublimationForceFactor*Mgas*grain_sublimation_rate_at_1AU*pow(r_h_AU,-2)*v_gas_h*grainArea;
        }
        
        if (1 && !gas_plot_run) {
            gmp_printf("%12.6f %12.3f %12.6f %12.6f %12.6f %12.6f %10.3f %10.3f %10.3f %10.6f\n",
                       orsa::FromUnits(t.get_d(),orsa::Unit::DAY,-1),
                       orsa::FromUnits(R_c.length(),orsa::Unit::KM,-1),
                       dist_ratio,
                       orsa::radToDeg()*acos(std::min(1.0,v_gas_g.normalized()*(R_c.normalized()))),
                       orsa::radToDeg()*acos(std::min(1.0,(vGrain-vComet).normalized()*(R_c.normalized()))),
                       orsa::radToDeg()*acos(std::min(1.0,(rSun-rComet).normalized()*(R_c.normalized()))),
                       orsa::FromUnits(R_c*uS,orsa::Unit::KM,-1),
                       orsa::FromUnits(R_c*uV,orsa::Unit::KM,-1),
                       orsa::FromUnits(R_c*uN,orsa::Unit::KM,-1),
                       orsa::FromUnits(grainRadius,orsa::Unit::METER,-1));
        }
        
        if (gas_plot_run) {
            // print for gas test
            const orsa::Vector r_local = g2l*(rGrain-rComet);
            const orsa::Vector v_gas_l = g2l*v_gas_g;
            const double test_grain_area = orsa::pi()*orsa::square(orsa::FromUnits(0.010000,orsa::Unit::METER));
            const orsa::Vector unit_thrust = 0.5*rho*v_gas_l.lengthSquared()*Cd*test_grain_area*v_gas_l.normalized();
            gmp_printf("GAS_TEST %+8.3f %+8.3f %+8.3f %8.3f   %8.3e   %+10.6f %+10.6f %+10.6f %10.6f   %+10.3e %+10.3e %+10.3e %10.3e\n",
                       orsa::FromUnits(r_local.getX(),orsa::Unit::KM,-1),
                       orsa::FromUnits(r_local.getY(),orsa::Unit::KM,-1),
                       orsa::FromUnits(r_local.getZ(),orsa::Unit::KM,-1),
                       orsa::FromUnits(r_local.length(),orsa::Unit::KM,-1),
                       //
                       number_density_gas_at_grain,
                       //
                       orsa::FromUnits(orsa::FromUnits(v_gas_l.getX(),orsa::Unit::KM,-1),orsa::Unit::SECOND),
                       orsa::FromUnits(orsa::FromUnits(v_gas_l.getY(),orsa::Unit::KM,-1),orsa::Unit::SECOND),
                       orsa::FromUnits(orsa::FromUnits(v_gas_l.getZ(),orsa::Unit::KM,-1),orsa::Unit::SECOND),
                       orsa::FromUnits(orsa::FromUnits(v_gas_l.length(),orsa::Unit::KM,-1),orsa::Unit::SECOND),
                       //
                       unit_thrust.getX(),
                       unit_thrust.getY(),
                       unit_thrust.getZ(),
                       unit_thrust.length());
            
        }
        
#warning make sure the expressions used are for a thrust, not for an acceleration
        
        // return thrust;
        return thrust+sublimationForce;
    }
public:
    bool nextEventTime(orsa::Time      &,
                       const mpz_class &) const {
        return false; // always ON
    }
protected:
    osg::ref_ptr<orsa::BodyGroup>  bg;
    osg::ref_ptr<const orsa::Body> sun;
    osg::ref_ptr<const orsa::Body> comet;
    osg::ref_ptr<const orsa::Body> grain;
    // const double grainBeta;
    const double grainDensity;
    // const double grainRadius;
    // const double grainArea;
    // const double Q_1AU; // gas production rate at 1 AU [units: number/second]
    const double eta_water; // nucleus_water_production_rate_factor
    const double Vgas_1AU; // gas velocity at 1 AU
    const double Mgas; // gas molecule mass (already converted from molar to KG)
    const double Cd; // drag coefficient
    const double newton;
};

class CGDIntegrator : public orsa::IntegratorRadau {
public:
    // gB: grain
    // nB: comet nucleus
    CGDIntegrator(const orsa::Body * gB,
                  const double & grain_initial_radius,
                  const double & grain_density,
                  const orsa::Body * nB,
                  const double & bound_distance,
                  const size_t & pow_10_max_distance,
                  const orsa::Time & t0_) :
        orsa::IntegratorRadau(),
        grain(gB),
        grainInitialRadius(grain_initial_radius),
        grainDensity(grain_density),
        nucleus(nB),
        r_bound(bound_distance),
        crossing_size(1+pow_10_max_distance),
        t0(t0_) {
        _accuracy = 1.0e-3;
        outcome = ORBITING;
        crossing_distance.resize(crossing_size);
        crossing_velocity.resize(crossing_size);
        crossing_time.resize(crossing_size);
        for (size_t k=0; k<crossing_size; ++k) {
            crossing_distance[k] = orsa::FromUnits(pow(10,k),orsa::Unit::KM);
            crossing_time[k] = t0;
            crossing_velocity[k] = -99.0;
        }
    }
protected:
    virtual ~CGDIntegrator() { }
public:
    enum OUTCOME_TYPE {
        ESCAPE=1,
        IMPACT=2,
        ORBITING=3
    };
protected:
    const orsa::Body * grain;
    const double grainInitialRadius;
    const double grainDensity;
    const orsa::Body * nucleus;
    const double r_bound;
public:
    mutable OUTCOME_TYPE outcome;
    mutable orsa::Cache<double> max_distance;
    const size_t crossing_size;
    const orsa::Time & t0;
    std::vector<double> crossing_distance;
    mutable std::vector<double> crossing_velocity;
    mutable std::vector<orsa::Time> crossing_time;
public:
    void singleStepDone(orsa::BodyGroup  * bg,
                        const orsa::Time & call_t,
                        const orsa::Time & call_dt,
                        orsa::Time       & ) const {
        const orsa::Time t = call_t + call_dt;
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
        max_distance.setIfLarger(grain_r_relative_local.length());
        orsa::IBPS nucleusIBPS;
        bg->getIBPS(nucleusIBPS,
                    nucleus, 
                    t);
        if (grain_r_relative_local.length() > r_bound) {
            // ORSA_DEBUG("escaped, aborting integration");
            outcome = ESCAPE;
            // ORSA_DEBUG("outcome: %i",outcome);
            // orsa::Integrator::abort();
        }
        
        for (size_t k=0; k<crossing_size; ++k) {
            if (grain_r_relative_local.length() > crossing_distance[k]) {
                if (crossing_time[k] == t0) {
                    crossing_time[k] = t;
                    crossing_velocity[k] = grain_v_relative_local.length();
                }
            }
        }
        
        if (nucleusIBPS.inertial->originalShape()->isInside(grain_r_relative_local)) {
            // ORSA_DEBUG("collision, aborting integration");
#warning note: the collision is not resolved exactly (i.e. rewind time for exact contact of body surface)
            outcome = IMPACT;
            // ORSA_DEBUG("outcome: %i",outcome);
            orsa::Integrator::abort();
        }
        
        if (grain_r_relative_local.length() > crossing_distance[crossing_size-1]) {
            // body reached largest tracked distance
            orsa::Integrator::abort();
        }
        
    }
};

#endif // COMET_GRAINS_DYNAMICS_H
