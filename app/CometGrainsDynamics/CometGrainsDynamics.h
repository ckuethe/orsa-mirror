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

using namespace orsa;

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
    GrainUpdateIBPS(const GrainUpdateIBPS & grain_ibps) {
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
                                     const double & sublimationRate, /* molecules per unit area per unit time */
                                     const double & moleculeMass, /* water molecule mass */
                                     const orsa::Body * grain) : 
        InertialBodyProperty(),
        _t0(t0),
        _initialRadius(initialRadius),
        _density(density),
        _dRadius_dt(sublimationRate*moleculeMass/(4.0*density)), /* check factors! i.e. the factor 4.0 due to energy balance, so it's sublimating from pi*Rg^2 instead of 4*pi*Rg^2 */
        _grain(grain) {
        // ORSA_DEBUG("dRdt: %g",_dRadius_dt);
        _init();
    }
public:
    GrainDynamicInertialBodyProperty(const GrainDynamicInertialBodyProperty & ibp) : 
        InertialBodyProperty(ibp),
        _t0(ibp._t0),
        _initialRadius(ibp._initialRadius),
        _density(ibp._density),
        _dRadius_dt(ibp._dRadius_dt),
        _grain(ibp._grain) {
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
    const double _dRadius_dt; // radius shrinking rate
    const orsa::Body * _grain;
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
public:
    InertialBodyProperty * clone() const {
        return new GrainDynamicInertialBodyProperty(*this);
    }
public:
    BodyPropertyType type() const { return BP_DYNAMIC; }
public:
    bool update(const orsa::Time & t) {
        
        // orsa::print(t);
        const double dt = (t-_t0).get_d();
        
        _radius = _initialRadius - _dRadius_dt * dt;
        
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
            const double & gas_production_rate_at_1AU,
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
        Q_1AU(gas_production_rate_at_1AU),
        Vgas_1AU(gas_velocity_at_1AU),
        Mgas(orsa::FromUnits(gas_molar_mass*1.66e-27,orsa::Unit::KG)), // conversion from molar
        Cd(gas_drag_coefficient),
        newton(orsa::FromUnits(orsa::FromUnits(orsa::FromUnits(1,orsa::Unit::KG),orsa::Unit::METER),orsa::Unit::SECOND,-2)) { }
protected:
    virtual ~GasDrag() { }
public:	
    orsa::Vector getThrust(const orsa::Time & t) const {
        
        // ORSA_DEBUG("referenceCount(): %i",referenceCount());
        
        if (Cd == 0.0) return orsa::Vector(0,0,0);
        
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
        
        const orsa::Vector R_h = (rComet-rSun);
        const double r_h = R_h.length();

        const double r_h_AU = orsa::FromUnits(r_h,orsa::Unit::AU,-1);
        
        const orsa::Vector R_c = (rGrain-rComet);
        const double r_c = R_c.length();
        
        // gas velocity at r_h, relative to comet
        const double v_gas_h = Vgas_1AU * pow(r_h_AU,-0.5);
        
        // production rate at r_h
        const double Q_h = Q_1AU * pow(r_h_AU,-2);
        
        // number density at r_c for production rate Q_h
        //const double n = Q_h / (4*orsa::pi()*orsa::square(r_c)*v_gas_h);
        // optional: can multiply x cos(theta_sun) to account for gas only from lit side of comet
        const double theta_sun = acos((rSun-rComet).normalized() * R_c.normalized());
        // const double theta_sun_factor = cos(0.5*theta_sun);
        const double theta_sun_factor = std::max(0.0,pow(cos(theta_sun),0.25)); // temperature for a low thermal inertia body goes as cos(theta_sun)^(1/4)
        const double n = Q_h * theta_sun_factor / (4*orsa::pi()*orsa::square(r_c)*v_gas_h);
        
        // ORSA_DEBUG("theta_sun_factor: %g",theta_sun_factor);
        
        // rho = mass density = number density x molecular mass
        const double rho = n * Mgas;
        
        const orsa::Matrix g2l = orsa::globalToLocal(comet,bg,t);
        const orsa::Matrix l2g = orsa::localToGlobal(comet,bg,t);
        orsa::IBPS ibps;
        if (!bg->getIBPS(ibps,comet,t)) {
            ORSA_DEBUG("problems at t:");
            orsa::print(t);
            orsa::crash();
        }
        const orsa::EllipsoidShape * nucleus_shape =
            dynamic_cast<const orsa::EllipsoidShape *> (ibps.inertial->originalShape());
        double na, nb, nc;
        nucleus_shape->getABC(na,nb,nc);
        const double nucleus_max_radius = std::max(na,std::max(nb,nc));
        
        if (nucleus_shape->isInside(g2l*(rGrain-rComet))) {
            // ORSA_DEBUG("zero thrust for grain inside comet nucleus...");
            return orsa::Vector(0,0,0);
        }
        
        // switch to radial when difference is approximately smaller than 1 deg
        orsa::Vector u_gas;
        /* #warning shoud replace this test with one more physically sound, comparing rotation period with gas expansion time to reach the distance
           if (r_c/nucleus_max_radius > 1000) {
           // radial
           u_gas = R_c.normalized();
           } else {
        */
        {
            // relative to comet
            // const orsa::Vector V_Gas_c   = v_gas_h * (rGrain-rComet).normalized();
            // modify V_gas_c to smoothly decrease near nucleus
            const orsa::Vector dr_g = rGrain-rComet;
            const orsa::Vector dr_l = g2l*dr_g;
            const orsa::Vector closest_point =
                nucleus_shape->closestVertex(dr_l);
            const orsa::Vector normal_l =
                nucleus_shape->normalVector(closest_point);
            const orsa::Vector normal_g = l2g*normal_l;
            u_gas = normal_g;
        }

        orsa::IBPS grain_ibps;
        // GrainIBPS grain_ibps;
        bg->getIBPS(grain_ibps,grain,t);
        osg::ref_ptr <GrainDynamicInertialBodyProperty> inertial = dynamic_cast <GrainDynamicInertialBodyProperty*> (grain_ibps.inertial.get());
        const double grainRadius = inertial->radius();
        const double grainArea = orsa::pi()*orsa::square(grainRadius);
        
        // NOTE: gas direction is proportional to normal_g, which gets close to radial at large distances
        
        const double dist_ratio = r_c / nucleus_max_radius;
        const double v_Gas_factor = dist_ratio/(1.0+dist_ratio); // goes from 0.5 near nucleus to 1.0 asymptotically
        const orsa::Vector V_Gas_c = v_Gas_factor * v_gas_h * u_gas;
        const orsa::Vector V_Grain_c = vGrain-vComet;
        const double dV = (V_Grain_c - V_Gas_c)*(V_Gas_c.normalized());
        const double sign = (dV>0) ? -1 : +1;
        const orsa::Vector thrust =
            sign*0.5*rho*dV*dV*Cd*grainArea*V_Gas_c.normalized();
        
        // force due to sublimation? acc = mol*Vgas*Z / (Rg*rho_grain)
        orsa::Vector sublimationForce(0,0,0);
        if (0) {
#warning MUST pass the parameters as arguments...
#warning what is the multiplicative factor in front?
            const double sublimationForceFactor = 0.01;
            const double grain_sublimation_rate = orsa::FromUnits(orsa::FromUnits(1.0e17,orsa::Unit::CM,-2),orsa::Unit::SECOND,-1);
            sublimationForce =
                // -uS*sublimationForceFactor*Mgas*grain_sublimation_rate*v_gas_h/(grainRadius*grainDensity); // this is just acc
                -uS*sublimationForceFactor*Mgas*grain_sublimation_rate*v_gas_h*grainArea;
        }
        
        if (0) {
            
            gmp_printf("%12.6f %12.3f %12.6f %12.6f %12.6f %12.6f %g %g %g %g\n",
                       orsa::FromUnits(t.get_d(),orsa::Unit::DAY,-1),
                       orsa::FromUnits(R_c.length(),orsa::Unit::KM,-1),
                       dist_ratio,
                       orsa::radToDeg()*acos(std::min(1.0,u_gas*(R_c.normalized()))),
                       orsa::radToDeg()*acos(std::min(1.0,(vGrain-vComet).normalized()*(R_c.normalized()))),
                       orsa::radToDeg()*acos(std::min(1.0,(rSun-rComet).normalized()*(R_c.normalized()))),
                       orsa::FromUnits(R_c*uS,orsa::Unit::KM,-1),
                       orsa::FromUnits(R_c*uV,orsa::Unit::KM,-1),
                       orsa::FromUnits(R_c*uN,orsa::Unit::KM,-1),
                       orsa::FromUnits(grainRadius,orsa::Unit::METER,-1));
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
    const double Q_1AU; // gas production rate at 1 AU [units: number/second]
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
