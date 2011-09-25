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

class GasDrag : public orsa::Propulsion {
public:
    GasDrag(orsa::BodyGroup  * bg_in,
            const orsa::Body * sun_in,
            const orsa::Body * comet_in,
            const orsa::Body * grain_in,
            const double & grain_beta_in,
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
        grainBeta(grain_beta_in),
        grainDensity(grain_density_in),
        grainRadius(GrainBetaToRadius(grainBeta,grainDensity)),
        grainArea(orsa::pi()*orsa::square(grainRadius)),
        Q_1AU(gas_production_rate_at_1AU),
        Vgas_1AU(gas_velocity_at_1AU),
        Mgas(orsa::FromUnits(gas_molar_mass*1.66e-27,orsa::Unit::KG)), // conversion from molar
        Cd(gas_drag_coefficient),
        newton(orsa::FromUnits(orsa::FromUnits(orsa::FromUnits(1,orsa::Unit::KG),orsa::Unit::METER),orsa::Unit::SECOND,-2)) { }
protected:
    ~GasDrag() { }
public:	
    orsa::Vector getThrust(const orsa::Time & t) const {
        
        if (Cd == 0.0) return orsa::Vector(0,0,0);
        
        orsa::Vector rSun;
        if (!bg->getInterpolatedPosition(rSun,sun.get(),t)) {
            ORSA_DEBUG("problems...");
        }	
        
        orsa::Vector rComet, vComet;
        if (!bg->getInterpolatedPosVel(rComet,vComet,comet.get(),t)) {
            ORSA_DEBUG("problems...");
        }
        
        orsa::Vector rGrain,vGrain;
        if (!bg->getInterpolatedPosVel(rGrain,vGrain,grain.get(),t)) {
            ORSA_DEBUG("problems...");
        }	
        
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
        const double theta_sun_factor = cos(0.5*theta_sun);
        const double n = Q_h * theta_sun_factor / (4*orsa::pi()*orsa::square(r_c)*v_gas_h);
        
        // ORSA_DEBUG("theta_sun_factor: %g",theta_sun_factor);
        
        // rho = mass density = number density x molecular mass
        const double rho = n * Mgas;
        
        const orsa::Matrix g2l = orsa::globalToLocal(comet,bg,t);
        const orsa::Matrix l2g = orsa::localToGlobal(comet,bg,t);
        orsa::IBPS ibps;
        bg->getIBPS(ibps,comet,t);
        const orsa::EllipsoidShape * nucleus_shape =
            dynamic_cast<const orsa::EllipsoidShape *> (ibps.inertial->originalShape());
        double na, nb, nc;
        nucleus_shape->getABC(na,nb,nc);
        const double nucleus_max_radius = std::max(na,std::max(nb,nc));
        
        // switch to radial when difference is approximately smaller than 1 deg
        orsa::Vector u_gas;
        if (r_c/nucleus_max_radius > 100) {
            // radial
            u_gas = R_c.normalized();
        } else {
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
        
        // NOTE: gas direction is proportional to normal_g, which gets close to radial at large distances
        
        const double dist_ratio = r_c / nucleus_max_radius;
        const double v_Gas_factor = dist_ratio/(1.0+dist_ratio); // goes from 0.5 near nucleus to 1.0 asymptotically
        const orsa::Vector V_Gas_c = v_Gas_factor * v_gas_h * u_gas;
        const orsa::Vector V_Grain_c = vGrain-vComet;
        const double dV = (V_Grain_c - V_Gas_c)*(V_Gas_c.normalized());
        const double sign = (dV>0) ? -1 : +1;
        const orsa::Vector thrust =
            sign*0.5*rho*dV*dV*Cd*grainArea*V_Gas_c.normalized();

        gmp_printf("%12.6f %12.3f %12.6f %12.6f %12.6f %12.6f %g %g %g\n",
                   orsa::FromUnits(t.get_d(),orsa::Unit::DAY,-1),
                   orsa::FromUnits(R_c.length(),orsa::Unit::KM,-1),
                   dist_ratio,
                   orsa::radToDeg()*acos(std::min(1.0,u_gas*(R_c.normalized()))),
                   orsa::radToDeg()*acos(std::min(1.0,(vGrain-vComet).normalized()*(R_c.normalized()))),
                   orsa::radToDeg()*acos(std::min(1.0,(rSun-rComet).normalized()*(R_c.normalized()))),
                   orsa::FromUnits(R_c.getX(),orsa::Unit::KM,-1),
                   orsa::FromUnits(R_c.getY(),orsa::Unit::KM,-1),
                   orsa::FromUnits(R_c.getZ(),orsa::Unit::KM,-1));
        
        return thrust;
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
    const double grainBeta;
    const double grainDensity;
    const double grainRadius;
    const double grainArea;
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
                  const orsa::Body * nB,
                  const double & bound_distance,
                  const size_t & pow_10_max_distance) :
        orsa::IntegratorRadau(),
        grain(gB),
        nucleus(nB),
        r_bound(bound_distance),
        crossing_size(1+pow_10_max_distance) {
        _accuracy = 1.0e-6;
        outcome = ORBITING;
        crossing_distance.resize(crossing_size);
        crossing_velocity.resize(crossing_size);
        crossing_time.resize(crossing_size);
        for (size_t k=0; k<crossing_size; ++k) {
            crossing_distance[k] = orsa::FromUnits(pow(10,k),orsa::Unit::KM);
            crossing_time[k] = orsa::Time(-99);
            crossing_velocity[k] = -99;
        }
    }
public:
    enum OUTCOME_TYPE {
        ESCAPE=1,
        IMPACT=2,
        ORBITING=3
    };
protected:
    const orsa::Body * grain;
    const orsa::Body * nucleus;
    const double r_bound;
public:
    mutable OUTCOME_TYPE outcome;
    mutable orsa::Cache<double> max_distance;
    const size_t crossing_size;
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
                if (crossing_time[k] == orsa::Time(-99)) {
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
