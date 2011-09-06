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
                  const orsa::Body * nB) :
        orsa::IntegratorRadau(),
        grain(gB),
        nucleus(nB) {
        _accuracy = 1.0e-6;
    }
protected:
    const orsa::Body * grain;
    const orsa::Body * nucleus;
public:
    void singleStepDone(orsa::BodyGroup  * bg,
                        const orsa::Time & call_t,
                        const orsa::Time & call_dt,
                        orsa::Time       & ) const {
        const orsa::Time t = call_t + call_dt;
        orsa::Vector r;
        bg->getInterpolatedPosition(r,
                                    nucleus,
                                    t);
        const orsa::Vector nucleus_r_global = r;
        bg->getInterpolatedPosition(r,
                                    grain,
                                    t);
        const orsa::Vector grain_r_relative_global = r - nucleus_r_global;
        const orsa::Matrix g2l = orsa::globalToLocal(nucleus,bg,t);
        const orsa::Vector grain_r_relative_local = g2l*grain_r_relative_global;
        orsa::IBPS nucleusIBPS;
        bg->getIBPS(nucleusIBPS,
                    nucleus, 
                    t);
        if (nucleusIBPS.inertial->originalShape()->isInside(grain_r_relative_local)) {
            ORSA_DEBUG("collision, aborting integration");
            orsa::Integrator::abort();
        }
        
    }
};

#endif // COMET_GRAINS_DYNAMICS_H
