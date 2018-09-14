#ifndef _VESTA_INTERIOR_ANALYTIC_H_
#define _VESTA_INTERIOR_ANALYTIC_H_

#include <orsa/util.h>
#include <orsa/statistic.h>

#include <orsaPDS/RadioScienceGravity.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_siman.h>

#include "CubicChebyshevMassDistribution.h"
#include "penalty.h"

#include "simplex.h"

#include "CCMD2SH.h"

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;

// GSL Simulated Annealing

/* how many points do we try before stepping */      
#define N_TRIES 100 // 100 // 200             

/* how many iterations for each T? */
#define ITERS_FIXED_T 100 // 200 // 100 // 1000 // 

/* max step size in random walk */
#define STEP_SIZE 1.0 // 50.0 // 1000.0

/* Boltzmann constant */
#define K 1.0                   

/* initial temperature */
// #define T_INITIAL 0.008     
#define T_INITIAL 0.001

/* damping factor for temperature */
#warning no damping??
#define MU_T 1.000 // 1.001 // 1.000 // 1.010 // 1.003      
#define T_MIN 1.0e-5 // 2.0e-6

gsl_siman_params_t params  = {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
                              K, T_INITIAL, MU_T, T_MIN};

class SIMAN_xp {
public:
    orsa::Cache<double> R0_plate;
    orsa::Cache<double> R0_gravity;
    orsa::Cache<double> bulkDensity;
    std::vector<orsa::Vector> rv;
    std::vector<double> hv;
    std::vector<double> dfb;
    std::vector<orsa::Vector> sv;
    orsa::Cache<double> min_rotation_period;
    orsa::Cache<double> max_rotation_period;
    orsa::Cache<double> max_axis_tilt;
    orsa::Cache<double> GM;
    orsa::Cache<size_t> SH_degree;
    orsa::Cache<size_t> T_degree;
    orsa::Cache<size_t> T_size;
    gsl_vector * cT0;
    gsl_vector * * uK;
    orsa::Cache<size_t> uK_size;
    std::vector<double> factor;
    orsa::Cache<double> sampling_factor;
    orsa::Cache<double> minimumDensity;
    orsa::Cache<double> maximumDensity;
    orsa::Cache<double> penaltyThreshold;
    osg::ref_ptr<const LayerData> layerData;
    osg::ref_ptr<orsa::Shape> shapeModel;
    osg::ref_ptr< SimplexIntegration<mpfr::mpreal> > si;
    orsa::Cache<orsa::Vector> sampled_CM;
};


void SIMAN_copy (void * source, void * dest) {
    SIMAN_xp * s = (SIMAN_xp *) source;
    SIMAN_xp * d = (SIMAN_xp *) dest;
    d->R0_plate            = s->R0_plate;
    d->R0_gravity          = s->R0_gravity;
    d->bulkDensity         = s->bulkDensity;
    d->rv                  = s->rv;
    d->hv                  = s->hv;
    d->dfb                 = s->dfb;
    d->sv                  = s->sv;
    d->min_rotation_period = s->min_rotation_period;
    d->max_rotation_period = s->max_rotation_period;
    d->max_axis_tilt       = s->max_axis_tilt;
    d->GM                  = s->GM;
    d->SH_degree           = s->SH_degree;
    d->T_degree            = s->T_degree;
    d->T_size              = s->T_size;
    d->cT0                 = s->cT0;
    d->uK                  = &(s->uK[0]);
    d->uK_size             = s->uK_size;
    d->factor              = s->factor;
    d->sampling_factor     = s->sampling_factor;
    d->minimumDensity      = s->minimumDensity;
    d->maximumDensity      = s->maximumDensity;
    d->penaltyThreshold    = s->penaltyThreshold;
    d->layerData           = s->layerData;
    d->shapeModel          = s->shapeModel;
    d->si                  = s->si;
    d->sampled_CM          = s->sampled_CM;
}

void * SIMAN_copy_construct (void * xp) {
    SIMAN_xp * d = new SIMAN_xp;
    SIMAN_copy(xp,d);
    return d;
}

void SIMAN_destroy (void * xp) {
    delete (SIMAN_xp *) xp;
}

static std::string CCMDF_output_filename;

std::string penalty_string_util(const std::string & type,
                                const double & val) {
    char line[4096];
    gmp_sprintf(line,"%s %+.6f",type.c_str(),val);
    return line;
}

// potential estimated at x; rv = mass sources position; dv = density of sources; vR2 = virtualRadius squared; vR3 = .. cube 
// OmegaSq = (2*pi/T)^2, T= rot. period
/*
double massClusterPotential(const orsa::Vector & x, const double & GM, const double & OmegaSq, const std::vector<orsa::Vector> & rv, const std::vector<double> & dv, const double & vR2, const double & vR3) {
    // all with positive sign convention (grav and rot potentials)
    double potential = 0;
    double sum_weight = 0;
    for (size_t k=0; k<rv.size(); ++k) {
        double weight = dv[k];
        double d2jk = (x-rv[k]).lengthSquared();
        // this version causes "varicella"-like bumps on geoid surface, better using good-old softening
        // if (d2jk>vR2) {
        //     potential += weight/sqrt(d2jk);
        // } else {
        //     potential += weight*(3*vR2-d2jk)/(2*vR3);
        // }
        // version with softening distance = vR
        potential += weight/sqrt(d2jk+vR2);
        sum_weight += weight;
    }
    if (sum_weight != 0) {
        potential /= sum_weight;
    }
    potential *= GM;
    // then the rotational component
    // const double Omega = orsa::twopi()/rotation_period;
    potential += 0.5*(orsa::square(x.getX())+orsa::square(x.getY()))*OmegaSq;
    //
    return potential;
}
*/

// potential estimated at r; rv = mass sources position; dv = density of sources; vR2 = virtualRadius squared; vR3 = .. cube 
// OmegaSq = (2*pi/T)^2, T= rot. period
double massClusterPotential(const orsa::Vector & r, const double & GM, const orsa::Vector & omega, const std::vector<orsa::Vector> & rv, const std::vector<double> & dv, const double & vR2, const double & vR3) {
    // all with positive sign convention (grav and rot potentials)
    double potential = 0;
    double sum_weight = 0;
    for (size_t k=0; k<rv.size(); ++k) {
        double weight = dv[k];
        double d2jk = (r-rv[k]).lengthSquared();
        // this version causes "varicella"-like bumps on geoid surface, better using good-old softening
        // if (d2jk>vR2) {
        //     potential += weight/sqrt(d2jk);
        // } else {
        //     potential += weight*(3*vR2-d2jk)/(2*vR3);
        // }
        // version with softening distance = vR
        potential += weight/sqrt(d2jk+vR2);
        sum_weight += weight;
    }
    if (sum_weight != 0) {
        potential /= sum_weight;
    }
    potential *= GM;
    // then the rotational component
    potential += 0.5*orsa::externalProduct(omega,r).lengthSquared();
    //
    return potential;
}

orsa::Vector massClusterAcceleration(const orsa::Vector & r, const double & GM, const orsa::Vector & omega, const std::vector<orsa::Vector> & rv, const std::vector<double> & dv, const double & vR2, const double & vR3, const double & ds) {
    const double p0 = massClusterPotential(r,                     GM,omega,rv,dv,vR2,vR3);
    const double pX = massClusterPotential(r+orsa::Vector(ds,0,0),GM,omega,rv,dv,vR2,vR3);
    const double pY = massClusterPotential(r+orsa::Vector(0,ds,0),GM,omega,rv,dv,vR2,vR3);
    const double pZ = massClusterPotential(r+orsa::Vector(0,0,ds),GM,omega,rv,dv,vR2,vR3);
    const double aX = (pX-p0)/ds;
    const double aY = (pY-p0)/ds;
    const double aZ = (pZ-p0)/ds;
    return orsa::Vector(aX,aY,aZ);
}

double E1(void * xp) {
    
    const bool verbose = false;
    
    SIMAN_xp * x = (SIMAN_xp *) xp;
    
    /* for (size_t q=0; q<x->factor.size(); ++q) {
       ORSA_DEBUG("factor[%02i] = %+12.6e",q,x->factor[q]);
       }
    */
    
    gsl_vector * cT = gsl_vector_alloc(x->T_size);
    gsl_vector_memcpy(cT,x->cT0);
    for (size_t b=0; b<x->uK_size; ++b) {
        for (size_t j=0; j<x->T_size; ++j) {
            gsl_vector_set(cT,j,gsl_vector_get(cT,j)+x->factor[b]*gsl_vector_get(x->uK[b],j));
        }
    }
    
    CubicChebyshevMassDistribution::CoefficientType coeff;
    CubicChebyshevMassDistribution::resize(coeff,x->T_degree); 

    for (unsigned int i=0; i<=x->T_degree; ++i) {
        for (unsigned int j=0; j<=x->T_degree-i; ++j) {
            for (unsigned int k=0; k<=x->T_degree-i-j; ++k) {
                if (i+j+k<=x->T_degree) {
                    const size_t index = CubicChebyshevMassDistribution::index(i,j,k);
                    coeff[i][j][k] = gsl_vector_get(cT,index);
                }
            }            
        }
    }
    /* osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
       new CubicChebyshevMassDistribution(coeff,x->bulkDensity,x->R0_plate,x->layerData);
    */
    osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
        new CubicChebyshevMassDistribution(coeff,x->R0_plate,x->layerData);
    
    osg::ref_ptr< orsa::Statistic<double> > stat = new orsa::Statistic<double>;
    orsa::Cache<double> minDensity, maxDensity;
    std::vector<double> dv;
    dv.resize(x->rv.size());
    for (size_t k=0; k<x->rv.size(); ++k) {
        dv[k] = massDistribution->density(x->rv[k]);
        stat->insert(dv[k]);
        minDensity.setIfSmaller(dv[k]);
        maxDensity.setIfLarger(dv[k]);
    }
    // const double minDensity = stat->min();
    // const double maxDensity = stat->max();
    const double averageSampledDensity = stat->average(); // can differ a bit from nominal average density
    
    // sample rotation values
    const double rotation_period = x->min_rotation_period + (x->max_rotation_period - x->min_rotation_period)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
    
    // rotation axis
    const double tilt = orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform() * x->max_axis_tilt;
    const double azim = (tilt != 0.0) ? orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform() * orsa::twopi() : 0.0;
    const orsa::Vector omega = (orsa::twopi()/rotation_period)*orsa::Vector(sin(tilt)*cos(azim),sin(tilt)*sin(azim),cos(tilt));
    
    // make sure this is called before leaving...
    gsl_vector_free(cT);
    
    std::vector<std::string> pv;
    pv.reserve(10);
    
    // start to fill pv with rotation values
    pv.push_back(penalty_string_util("TROT",orsa::FromUnits(rotation_period,orsa::Unit::HOUR,-1)));
    pv.push_back(penalty_string_util("TILT",orsa::radToDeg()*tilt));
    pv.push_back(penalty_string_util("AZIM",orsa::radToDeg()*azim));
    
    // condensed all variants here below...
    /*
    double retVal =
        1000.0*std::max(0.0,(x->minimumDensity-minDensity)) +
        1000.0*std::max(0.0,(maxDensity-x->maximumDensity));
    */
    double retVal = 0.0;
    if ((minDensity<x->minimumDensity) || (maxDensity>x->maximumDensity)) {
        retVal += 1e6;
        if (minDensity<x->minimumDensity) {
            retVal += 1e6*(x->minimumDensity-minDensity);
        }
        if (maxDensity>x->maximumDensity) {
            retVal += 1e6*(maxDensity-x->maximumDensity);
        }
    }
    double penalty = retVal;
    pv.push_back(penalty_string_util("DR",penalty));
    if ( (minDensity > x->minimumDensity) && (maxDensity < x->maximumDensity) ) {
        
        // target: low penalty value...
        
        if (0) {
            // target: no low-density "holes"
            double delta_penalty = 0.0;
            size_t entries=0;
            for (size_t k=1; k<dv.size(); ++k) {
                // assuming the mid-point rm is still contained in body...
                const orsa::Vector rm = 0.5*(x->rv[k] + x->rv[k-1]);
                const double dm = massDistribution->density(rm);
                const double d12 = std::min(dv[k],dv[k-1]);
                if (dm<d12) {
                    delta_penalty += (d12-dm)/x->bulkDensity * x->R0_plate/(x->rv[k] - x->rv[k-1]).length();
                    ++entries;
                }
            }
            if (entries!=0) delta_penalty /= entries;
            delta_penalty *= 1000.0;
            
            penalty += delta_penalty;
            // pv.push_back(penalty_string_util("NDH",delta_penalty));
            pv.push_back(penalty_string_util("CDP",delta_penalty)); // Convex Density Profile
            if (verbose) ORSA_DEBUG("delta penalty: %+10.6f   [target: no low-density \"holes\"]",delta_penalty);
        }
                
        if (1 && (x->sv.size()>0)) {
            // XP: for points of given density, test value of potential
            const double virtualRadius = cbrt(x->shapeModel->volume()/(x->rv.size()*(4.0*orsa::pi()/3.0)));
            const double virtualRadiusSq = virtualRadius*virtualRadius;
            const double virtualRadiusCb = virtualRadius*virtualRadiusSq;
            // const double OmegaSq = orsa::square(orsa::twopi()/rotation_period);
            const double km   = orsa::FromUnits(1.0,orsa::Unit::KM);
            const double gcm3 = orsa::FromUnits(orsa::FromUnits(1.0,orsa::Unit::GRAM),orsa::Unit::CM,-3);
            const double delta_rr = 1000.0*km; // 50.0*km; // 100.0*km;
            double rr = 0.0*km; // ref. radius
            while (1) {
                rr += delta_rr;
                osg::ref_ptr< orsa::Statistic<double> > stat_rho = new orsa::Statistic<double>;
                for (size_t j=0; j<x->sv.size(); ++j) {
                    const orsa::Vector & v = x->sv[j];
                    if (orsa::square(rr) < v.lengthSquared()) {
                        const orsa::Vector u = v.normalized();
                        stat_rho->insert(massDistribution->density(u*rr));
                    }
                }
                if (stat_rho->entries() != x->sv.size()) {
                    // some test points were outside the shape...
                    break;
                }
                const double rhoRef = stat_rho->average();
                //
                osg::ref_ptr< orsa::Statistic<double> > stat_xp = new orsa::Statistic<double>;
                // std::vector<double> xp; // potential values
                osg::ref_ptr< orsa::Statistic<double> > stat_ll = new orsa::Statistic<double>;
                for (size_t j=0; j<x->sv.size(); ++j) {
                    // search (bisection) for point with density rhoRef along sv[j];
                    const orsa::Vector u = x->sv[j].normalized();
                    double lA = 0.0;
                    double rhoA = massDistribution->density(u*lA);
                    double lB = x->sv[j].length();
                    double rhoB = massDistribution->density(u*lB);
                    bool converged=false;
                    if ((rhoRef>=std::min(rhoA,rhoB)) && (rhoRef<=std::max(rhoA,rhoB))) {
                        while (1) {
                            // ORSA_DEBUG("lA: %g rhoA: %g   lB: %g rhoB: %g",lA/km,rhoA/gcm3,lB/km,rhoB/gcm3);
                            if (fabs(lA-lB)<orsa::FromUnits(0.1,orsa::Unit::KM)) { converged=true; break; }
                            double lC = 0.5*(lA+lB);
                            double rhoC = massDistribution->density(u*lC);
                            if ((rhoRef>=std::min(rhoA,rhoC)) && (rhoRef<=std::max(rhoA,rhoC))) {
                                lB   = lC;
                                rhoB = rhoC;
                            } else {
                                lA   = lC;
                                rhoA = rhoC;
                            }
                        }
                    }
                    if (converged) {
                        double lC = 0.5*(lA+lB);
                        double potential = massClusterPotential(u*lC,x->GM,omega,x->rv,dv,virtualRadiusSq,virtualRadiusCb);
                        stat_xp->insert(potential);
                        // xp.push_back(potential);
                        stat_ll->insert(lC);
                    }
                }
                // double delta_penalty = stat_xp->variance();
                // delta_penalty *= 0.0001; // norm.
                //
                double delta_penalty = stat_xp->standardDeviation()/stat_xp->average();
                if (!finite(delta_penalty)) {
                    delta_penalty = 1.0;
                }
                delta_penalty *= 10000.0; // norm.
                //
                penalty += delta_penalty;
                // ORSA_DEBUG("XP vR: %g   V: %g +/- %g",virtualRadius,stat_xp->average(),stat_xp->averageError());
                char name[4096];
                // sprintf(name,"XP%.1f@%g",rhoRef/gcm3,rr/km);
                sprintf(name,"XP%.2f@%g:%.1f:%.1f",rhoRef/gcm3,rr/km,stat_ll->min()/km,stat_ll->max()/km);
                pv.push_back(penalty_string_util(name,delta_penalty));
            }
        }
        
        if (0 && (x->sv.size()>0)) {
            // RP: range of potential, minimum and maximum radius of points having the same potential
            // output only, not included in penalty computation, not need in production...
            const double virtualRadius = cbrt(x->shapeModel->volume()/(x->rv.size()*(4.0*orsa::pi()/3.0)));
            const double virtualRadiusSq = virtualRadius*virtualRadius;
            const double virtualRadiusCb = virtualRadius*virtualRadiusSq;
            // const double OmegaSq = orsa::square(orsa::twopi()/rotation_period);
            const double km   = orsa::FromUnits(1.0,orsa::Unit::KM);
            const double gcm3 = orsa::FromUnits(orsa::FromUnits(1.0,orsa::Unit::GRAM),orsa::Unit::CM,-3);
            const double delta_rr = 100.0*km; // 50.0*km; // 100.0*km;
            double rr = 0.0*km; // ref. radius
            while (1) {
                rr += delta_rr;
                osg::ref_ptr< orsa::Statistic<double> > stat_potential = new orsa::Statistic<double>;
                for (size_t j=0; j<x->sv.size(); ++j) {
                    const orsa::Vector & v = x->sv[j];
                    if (orsa::square(rr) < v.lengthSquared()) {
                        const orsa::Vector u = v.normalized();
                        stat_potential->insert(massClusterPotential(u*rr,x->GM,omega,x->rv,dv,virtualRadiusSq,virtualRadiusCb));
                    }
                }
                if (stat_potential->entries() != x->sv.size()) {
                    // some test points were outside the shape...
                    break;
                }
                const double potentialRef = stat_potential->average();
                //
                // osg::ref_ptr< orsa::Statistic<double> > stat_rp = new orsa::Statistic<double>;
                // std::vector<double> rp; // potential values
                osg::ref_ptr< orsa::Statistic<double> > stat_ll = new orsa::Statistic<double>;
                for (size_t j=0; j<x->sv.size(); ++j) {
                    // search (bisection) for point with potentialRef along sv[j];
                    const orsa::Vector u = x->sv[j].normalized();
                    double lA = 0.0;
                    double potentialA = massClusterPotential(u*lA,x->GM,omega,x->rv,dv,virtualRadiusSq,virtualRadiusCb);
                    double lB = x->sv[j].length();
                    double potentialB = massClusterPotential(u*lB,x->GM,omega,x->rv,dv,virtualRadiusSq,virtualRadiusCb);
                    bool converged=false;
                    if ((potentialRef>=std::min(potentialA,potentialB)) && (potentialRef<=std::max(potentialA,potentialB))) {
                        while (1) {
                            // ORSA_DEBUG("lA: %g potentialA: %g   lB: %g potentialB: %g",lA/km,potentialA/gcm3,lB/km,potentialB/gcm3);
                            if (fabs(lA-lB)<orsa::FromUnits(0.1,orsa::Unit::KM)) { converged=true; break; }
                            double lC = 0.5*(lA+lB);
                            double potentialC = massClusterPotential(u*lC,x->GM,omega,x->rv,dv,virtualRadiusSq,virtualRadiusCb);
                            if ((potentialRef>=std::min(potentialA,potentialC)) && (potentialRef<=std::max(potentialA,potentialC))) {
                                lB   = lC;
                                potentialB = potentialC;
                            } else {
                                lA   = lC;
                                potentialA = potentialC;
                            }
                        }
                    }
                    if (converged) {
                        double lC = 0.5*(lA+lB);
                        double potential = massClusterPotential(u*lC,x->GM,omega,x->rv,dv,virtualRadiusSq,virtualRadiusCb);
                        // stat_rp->insert(potential);
                        // rp.push_back(potential);
                        stat_ll->insert(lC);
                    }
                }
                // double delta_penalty = stat_rp->variance();
                // delta_penalty *= 0.0001; // norm.
                
                // double delta_penalty = stat_rp->standardDeviation()/stat_rp->average();
                // delta_penalty *= 10000.0; // norm.
                //
                // penalty += delta_penalty;
                // ORSA_DEBUG("RP vR: %g   V: %g +/- %g",virtualRadius,stat_rp->average(),stat_rp->averageError());
                const double delta_penalty = 0.0;
                
                char name[4096];
                // sprintf(name,"RP%.1f@%g",rhoRef/gcm3,rr/km);
                // sprintf(name,"RP%.1f@%g:%.1f:%.1f",potentialRef,rr/km,stat_ll->min()/km,stat_ll->max()/km);
                // sprintf(name,"RP%.1f@%g:%.1f:%.1f",potentialRef,rr/km,stat_ll->min()/km,stat_ll->max()/km);
                sprintf(name,"RP@%g:%.1f:%.1f",rr/km,stat_ll->min()/km,stat_ll->max()/km);
                pv.push_back(penalty_string_util(name,delta_penalty));
            }
        }
        
        if (0 && (x->sv.size()>0)) {
            // SP: test value of potential at surface
            const double virtualRadius = cbrt(x->shapeModel->volume()/(x->rv.size()*(4.0*orsa::pi()/3.0)));
            const double virtualRadiusSq = virtualRadius*virtualRadius;
            const double virtualRadiusCb = virtualRadius*virtualRadiusSq;
            // const double OmegaSq = orsa::square(orsa::twopi()/rotation_period);
            const double km   = orsa::FromUnits(1.0,orsa::Unit::KM);
            const double gcm3 = orsa::FromUnits(orsa::FromUnits(1.0,orsa::Unit::GRAM),orsa::Unit::CM,-3);
            osg::ref_ptr< orsa::Statistic<double> > stat_sp = new orsa::Statistic<double>;
            // osg::ref_ptr< orsa::Statistic<double> > stat_ll = new orsa::Statistic<double>;
            for (size_t j=0; j<x->sv.size(); ++j) {
                double potential = massClusterPotential(x->sv[j],x->GM,omega,x->rv,dv,virtualRadiusSq,virtualRadiusCb);
                stat_sp->insert(potential);
                // stat_ll->insert(x->sv[j].length());
            }
            //
            double delta_penalty = stat_sp->standardDeviation()/stat_sp->average();
            delta_penalty *= 10000.0; // norm.
            //
            penalty += delta_penalty;
            // ORSA_DEBUG("SP vR: %g   V: %g +/- %g",virtualRadius,stat_sp->average(),stat_sp->averageError());
            char name[4096];
            sprintf(name,"SP");
            pv.push_back(penalty_string_util(name,delta_penalty));
        }
        
        if (0) {
            // target: density proportional to depth
            double delta_penalty = 0.0;
            for (size_t k=0; k<dv.size(); ++k) {
                delta_penalty += (1.0 - dv[k]/x->bulkDensity)*(1.0 + x->hv[k]/x->R0_plate);
            }
            delta_penalty /= dv.size();
            delta_penalty *= 1.0;
            penalty += delta_penalty;
            pv.push_back(penalty_string_util("DID",delta_penalty));
            if (verbose) ORSA_DEBUG("delta penalty: %+10.6f   [target: density proportional to depth]",delta_penalty);
        }
        
        if (0) {
            // target: density decreasing from barycenter
            double delta_penalty = 0.0;
            for (size_t k=0; k<dv.size(); ++k) {
                delta_penalty += (1.0 - dv[k]/x->bulkDensity)*(1.0 - x->dfb[k]/x->R0_plate);
            }
            delta_penalty /= dv.size();
            delta_penalty *= 100.0;
            penalty += delta_penalty;
            pv.push_back(penalty_string_util("DDB",delta_penalty));
            if (verbose) ORSA_DEBUG("delta penalty: %+10.6f   [target: density decreasing from barycenter]",delta_penalty);
        }
                
        if (0) {
            // target: no low-density "holes", with a SCALE
            const double max_distance = orsa::FromUnits(5.0,orsa::Unit::KM);
            const double max_d2 = max_distance*max_distance;
            double delta_penalty = 0.0;
            size_t entries=0;
            for (size_t k=1; k<dv.size(); ++k) {
                if ((x->rv[k-1]-x->rv[k]).lengthSquared() > max_d2) continue;
                // assuming the mid-point rm is still contained in body...
                const orsa::Vector rm = 0.5*(x->rv[k] + x->rv[k-1]);
                const double dm = massDistribution->density(rm);
                const double d12 = std::min(dv[k],dv[k-1]);
                if (dm<d12) {
                    delta_penalty += (d12-dm)/x->bulkDensity * x->R0_plate/(x->rv[k] - x->rv[k-1]).length();
                    ++entries;
                }
            }
            delta_penalty *= 200.0;
            if (entries!=0) delta_penalty /= entries;
            
            penalty += delta_penalty;
            // pv.push_back(penalty_string_util("NDH",delta_penalty));
            pv.push_back(penalty_string_util("CDPS",delta_penalty)); // Convex Density Profile
            if (verbose) ORSA_DEBUG("delta penalty: %+10.6f   [target: no low-density \"holes\"]",delta_penalty);
        }
        
        if (0) {
            // target: closest to uniform (simple)
            const double delta_penalty = (maxDensity-minDensity)/x->bulkDensity;
            penalty += delta_penalty;
            pv.push_back(penalty_string_util("MINDR",delta_penalty));
            if (verbose) ORSA_DEBUG("delta penalty: %+10.6f   [target: closest to uniform (simple)]",delta_penalty);
        }
        
        if (0) {
            // target: highest single density peak
            const double delta_penalty = (minDensity-maxDensity)/x->bulkDensity;
            penalty += delta_penalty;
            pv.push_back(penalty_string_util("MAXDR",delta_penalty));
            if (verbose) ORSA_DEBUG("delta penalty: %+10.6f   [target: highest single density peak]",delta_penalty);
        }
        
        if (0) {
            // target: maximum moment of inertia Izz = C/MR^2
            orsa::Cache<orsa::Vector> CM = x->sampled_CM;
            mpf_class IxxMR2, IyyMR2, IzzMR2;
            inertia(CM,IxxMR2,IyyMR2,IzzMR2,x->si.get(),massDistribution.get(),x->R0_plate);            
            const double delta_penalty = -5000.0*IzzMR2.get_d(); // -5000
            penalty += delta_penalty;
            pv.push_back(penalty_string_util("MAXIZZ",delta_penalty));
        }
        
        if (0) {
            // target: minimum moment of inertia Izz = C/MR^2
            orsa::Cache<orsa::Vector> CM = x->sampled_CM;
            mpf_class IxxMR2, IyyMR2, IzzMR2;
            inertia(CM,IxxMR2,IyyMR2,IzzMR2,x->si.get(),massDistribution.get(),x->R0_plate);            
            const double delta_penalty = +5000.0*IzzMR2.get_d(); // +5000
            penalty += delta_penalty;
            pv.push_back(penalty_string_util("MINIZZ",delta_penalty));
        }
        
        if (0) {
            // surface density standard deviation (pushes towards uniform value of surface density)
            osg::ref_ptr< orsa::Statistic<double> > stat_sd = new orsa::Statistic<double>;
            for (size_t k=0; k<x->sv.size(); ++k) {
                stat_sd->insert(massDistribution->density(x->sv[k]));
            }
            double delta_penalty = stat_sd->standardDeviation()/x->bulkDensity;
            delta_penalty *= 10.0;
            penalty += delta_penalty;
            pv.push_back(penalty_string_util("SD",delta_penalty));            
        }
        
        retVal = penalty;
    }
    
#warning do not call return here, call it only at the very end of this method

    /* 
    // OLD VERSIONS...
    // most flat
    // return (maxDensity-minDensity)+10000*(penalty/x->penaltyThreshold)+10*std::max(0.0,(x->minimumDensity-minDensity))+10*std::max(0.0,(maxDensity-x->maximumDensity));
    
    // most peaks
    // return (minDensity-maxDensity)+10000*(penalty/x->penaltyThreshold)+10*std::max(0.0,(x->minimumDensity-minDensity))+10*std::max(0.0,(maxDensity-x->maximumDensity));
    
    // generic
    // return 10000*(penalty/x->penaltyThreshold)+10*std::max(0.0,(x->minimumDensity-minDensity))+10*std::max(0.0,(maxDensity-x->maximumDensity));
    
    // most flat, no penalty
    // return (maxDensity-minDensity)+10*std::max(0.0,(x->minimumDensity-minDensity))+10*std::max(0.0,(maxDensity-x->maximumDensity));
    
    // most peaks, no penalty
    // return (minDensity-maxDensity)+10*std::max(0.0,(x->minimumDensity-minDensity))+10*std::max(0.0,(maxDensity-x->maximumDensity));
    
    // generic, no penalty
    // return 10*std::max(0.0,(x->minimumDensity-minDensity))+10*std::max(0.0,(maxDensity-x->maximumDensity));
    
    // most peaks, negative penalty
    // return (minDensity-maxDensity)-10000*(penalty/x->penaltyThreshold)+10*std::max(0.0,(x->minimumDensity-minDensity))+10*std::max(0.0,(maxDensity-x->maximumDensity));
    */
    
    char pvline[4096];
    sprintf(pvline,"");
    {
        char str[4096];
        for (size_t p=0; p<pv.size(); ++p) {
            if (p==0) {
                sprintf(str, "%s",pv[p].c_str());
            } else {
                sprintf(str," %s",pv[p].c_str());
            }
            strcat(pvline,str);
        }
    }
    // ORSA_DEBUG("pvline: [%s]",pvline);
    
    ORSA_DEBUG("[density] min: %+6.2f max: %+6.2f avg: %+6.2f [g/cm^3]   penalty: %+10.6f [%s]",
               orsa::FromUnits(orsa::FromUnits(minDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
               orsa::FromUnits(orsa::FromUnits(maxDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
               orsa::FromUnits(orsa::FromUnits(averageSampledDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
               penalty,
               pvline);
    
    if ( (minDensity >= x->minimumDensity) &&
         (maxDensity <= x->maximumDensity) &&
         (penalty <= x->penaltyThreshold) ) {
        // another quick output...
#warning pass filename as parameter...
        CubicChebyshevMassDistributionFile::CCMDF_data data;
        data.minDensity = minDensity;
        data.maxDensity = maxDensity;
        data.deltaDensity = maxDensity-minDensity;
        data.penalty = penalty;
        // data.densityScale = x->bulkDensity;
        data.R0 = x->R0_plate;
        data.SH_degree = x->SH_degree;
        data.coeff = coeff;
        data.layerData = x->layerData;
        //
        if (1) {
            orsa::Cache<orsa::Vector> CM = x->sampled_CM;
            mpf_class IxxMR2, IyyMR2, IzzMR2;
            inertia(CM,
                    IxxMR2,
                    IyyMR2,
                    IzzMR2,
                    // x->SH_degree,
                    x->si.get(),
                    massDistribution.get(),
                    x->R0_plate);
            
            char comment[4096];
            sprintf(comment,"%.6f %.6f %.6f",IxxMR2.get_d(),IyyMR2.get_d(),IzzMR2.get_d());
            if (0) {
                char tmpstr[4096];
                for (size_t b=0; b<x->uK_size; ++b) {
                    sprintf(tmpstr," %+12.6f",x->factor[b]);
                    strcat(comment,tmpstr);
                }
            }
            strcat(comment," ");
            strcat(comment,pvline);
            data.comment = comment;
        }
        //
        // CubicChebyshevMassDistributionFile::append(data,"CCMDF.out");
        CubicChebyshevMassDistributionFile::append(data,CCMDF_output_filename.c_str());
    }
    
    return retVal;
}

double M1(void * xp, void * yp) {
    SIMAN_xp * x = (SIMAN_xp *) xp;
    SIMAN_xp * y = (SIMAN_xp *) yp;

    double distance = 0.0;
    for (size_t b=0; b<x->uK_size; ++b) {
        distance += orsa::square(y->factor[b]-x->factor[b]);
    }
    distance /= x->uK_size;
    distance = sqrt(distance);

    return distance;
}

void S1(const gsl_rng * r, void * xp, double step_size) {
    SIMAN_xp * x = (SIMAN_xp *) xp;
    
    for (size_t b=0; b<x->uK_size; ++b) {
        // x->factor[b] += step_size*(2*gsl_rng_uniform(r)-1)/x->uK_size;
        // x->factor[b] += step_size*(2*gsl_rng_uniform(r)-1);
        x->factor[b] += x->sampling_factor*step_size*(2*gsl_rng_uniform(r)-1);
    }
}

void P1(void *) {
    ORSA_DEBUG("print here...");
}



#endif // _VESTA_INTERIOR_ANALYTIC_H_
