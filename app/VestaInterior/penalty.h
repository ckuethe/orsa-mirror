#ifndef ORSA_INTERIOR_PENALTY_H
#define ORSA_INTERIOR_PENALTY_H

#include <orsa/util.h>

double MassDistributionPenalty(const std::vector<orsa::Vector> & rv,
                               const std::vector<double> & dv,
                               const orsa::MassDistribution * md) {
    double penalty = 0.0;
    for (size_t k1=0; k1<rv.size(); ++k1) {
        for (size_t k2=0; k2<k1; ++k2) {
#warning assuming all middle points are inside the body; slightly more complicated algorithm needed with stongly concave bodies
            const orsa::Vector rm = 0.5*(rv[k1]+rv[k2]);
            const double dm = md->density(rm);
            const double d12 = std::min(dv[k1],dv[k2]);
            if (dm<d12) {
                // penalty += (d12-dm)/d12;
                penalty = std::max(penalty,(d12-dm)/d12);
            }
        }
    }
    return penalty;
}

double MassDistributionPenalty(const std::vector<orsa::Vector> & rv,
                               const orsa::MassDistribution * md) {
    std::vector<double> dv;
    dv.resize(rv.size());
    for (size_t k=0; k<rv.size(); ++k) {
        dv[k] = md->density(rv[k]);
    }
    return MassDistributionPenalty(rv,dv,md);
}

double MassDistributionPenalty(const orsa::RandomPointsInShape * randomPointsInShape) {
    std::vector<orsa::Vector> rv;
    std::vector<double> dv;
    orsa::Vector v;
    double density;
    randomPointsInShape->reset();
    while (randomPointsInShape->get(v,density)) { 
        rv.push_back(v);
        dv.push_back(density);
    }
    return MassDistributionPenalty(rv,dv,randomPointsInShape->md.get());
}

#endif // ORSA_INTERIOR_PENALTY_H
