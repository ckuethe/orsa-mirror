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
        _accuracy = 1.0e-3;
        outcome = ORBITING;
        crossing_distance.resize(crossing_size);
        crossing_time.resize(crossing_size);
        for (size_t k=0; k<crossing_size; ++k) {
            crossing_distance[k] = orsa::FromUnits(exp10(k),orsa::Unit::KM);
            crossing_time[k] = orsa::Time(0);
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
                if (crossing_time[k] == orsa::Time(0)) {
                    crossing_time[k] = t;
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
