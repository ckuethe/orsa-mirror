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

/*** CHOOSE ONE ***/
// typedef double simplex_T;
// typedef mpf_class simplex_T;
typedef dd_real simplex_T;
// typedef qd_real simplex_T;

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;

// GSL Simulated Annealing

/* how many points do we try before stepping */      
#define N_TRIES 100 // 100 // 200             

/* how many iterations for each T? */
#define ITERS_FIXED_T 100 // 200 // 100 // 1000 // 

/* max step size in random walk */
#define STEP_SIZE 1000.0

/* Boltzmann constant */
#define K 1.0                   

/* initial temperature */
// #define T_INITIAL 0.008     
#define T_INITIAL 0.001

/* damping factor for temperature */
#warning no damping??
#define MU_T 1.001 // 1.000 // 1.010 // 1.003      
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
    orsa::Cache<size_t> SH_degree;
    orsa::Cache<size_t> T_degree;
    orsa::Cache<size_t> T_size;
    gsl_vector * cT0;
    gsl_vector * * uK;
    orsa::Cache<size_t> uK_size;
    std::vector<double> factor;
    orsa::Cache<double> minimumDensity;
    orsa::Cache<double> maximumDensity;
    orsa::Cache<double> penaltyThreshold;
    osg::ref_ptr<const LayerData> layerData;
    osg::ref_ptr<orsa::Shape> shapeModel;
    osg::ref_ptr< SimplexIntegration<simplex_T> > si;
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
    d->SH_degree           = s->SH_degree;
    d->T_degree            = s->T_degree;
    d->T_size              = s->T_size;
    d->cT0                 = s->cT0;
    d->uK                  = &(s->uK[0]);
    d->uK_size             = s->uK_size;
    d->factor              = s->factor;
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

std::string penalty_string_util(const std::string & type,
                                const double & val) {
    char line[4096];
    gmp_sprintf(line,"%s %+.6f",type.c_str(),val);
    return line;
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
    
    // make sure this is called before leaving...
    gsl_vector_free(cT);
    
    std::vector<std::string> pv;
    pv.reserve(10);
    
    // condensed all variants here below...
    double retVal  =
        10.0*std::max(0.0,(x->minimumDensity-minDensity)) +
        10.0*std::max(0.0,(maxDensity-x->maximumDensity));
    double penalty = retVal;
    pv.push_back(penalty_string_util("DR",penalty));
    if ( (minDensity > x->minimumDensity) &&
         (maxDensity < x->maximumDensity) ) {
        
        // target: low penalty value...
        
        if (0) {
            double delta_penalty = 0.0;
            // target: closest to uniform density
            for (size_t k=0; k<dv.size(); ++k) {
                delta_penalty += fabs(dv[k]/x->bulkDensity - 1.0);
            }
            delta_penalty /= dv.size();
            delta_penalty *= 10.0;
            penalty += delta_penalty;
            pv.push_back(penalty_string_util("ALT_MINDR",delta_penalty));
            if (verbose) ORSA_DEBUG("delta penalty: %+10.6f   [target: closest to uniform density]",delta_penalty);
        }
        
        if (0) {
            // target: closest to uniform (simple)
            const double delta_penalty = (maxDensity-minDensity)/x->bulkDensity;
            penalty += delta_penalty;
            pv.push_back(penalty_string_util("MINDR",delta_penalty));
            if (verbose) ORSA_DEBUG("delta penalty: %+10.6f   [target: closest to uniform (simple)]",delta_penalty);
        }
        
        if (0) { /* ALT */
            // target: most volume with high density
            // const double exponent = 1.0;
            double delta_penalty = 0.0;
            for (size_t k=0; k<dv.size(); ++k) {
                // v1
                // delta_penalty -= pow(dv[k]/x->bulkDensity,exponent);
                delta_penalty += 1.0 - dv[k]/x->bulkDensity;
                
                // v2 (which does the opposite of what it should...)
                /* if ((dv[k] - x->bulkDensity) > 0.85*(maxDensity - x->bulkDensity)) {
                   delta_penalty -= 1.0;
                   }
                */
            }
            delta_penalty /= dv.size();
            delta_penalty *= 1.0e3/sqrt(dv.size());
            penalty += delta_penalty;
            pv.push_back(penalty_string_util("MVHD",delta_penalty));
            if (verbose) ORSA_DEBUG("delta penalty: %+10.6f   [target: most volume with high density]",delta_penalty);
        }
        
        if (0) {
            // target: highest single density peak
            const double delta_penalty = (minDensity-maxDensity)/x->bulkDensity;
            penalty += delta_penalty;
            pv.push_back(penalty_string_util("MAXDR",delta_penalty));
            if (verbose) ORSA_DEBUG("delta penalty: %+10.6f   [target: highest single density peak]",delta_penalty);
        }
        
        if (0) {
#warning need to look more into this
            // target: fraction of mass at given density ranges
        }
        
        if (0) {
            // target: density proportional to depth
            double delta_penalty = 0.0;
            for (size_t k=0; k<dv.size(); ++k) {
                delta_penalty += (1.0 - dv[k]/x->bulkDensity)*(1.0 + x->hv[k]/x->R0_plate);
            }
            delta_penalty /= dv.size();
            delta_penalty *= 100.0;
            penalty += delta_penalty;
            pv.push_back(penalty_string_util("DID",delta_penalty));
            if (verbose) ORSA_DEBUG("delta penalty: %+10.6f   [target: density proportional to depth]",delta_penalty);
        }
        
        if (1) {
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
        
        if (1) {
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
            delta_penalty *= 100.0;
            if (entries!=0) delta_penalty /= entries;
            
            penalty += delta_penalty;
            pv.push_back(penalty_string_util("NDH",delta_penalty));
            if (verbose) ORSA_DEBUG("delta penalty: %+10.6f   [target: no low-density \"holes\"]",delta_penalty);
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
            orsa::Cache<orsa::Vector> CM;
            mpf_class IxxMR2, IyyMR2, IzzMR2;
            inertia(CM,
                    IxxMR2,
                    IyyMR2,
                    IzzMR2,
                    x->SH_degree,
                    x->si.get(),
                    massDistribution.get(),
                    x->R0_plate);
            
            char comment[4096];
            sprintf(comment,"%.6f %.6f %.6f",IxxMR2.get_d(),IyyMR2.get_d(),IzzMR2.get_d());
            char tmpstr[4096];
            for (size_t b=0; b<x->uK_size; ++b) {
                sprintf(tmpstr," %+12.6f",x->factor[b]);
                strcat(comment,tmpstr);
            }
            strcat(comment," ");
            strcat(comment,pvline);
            data.comment = comment;
        }
        //
        CubicChebyshevMassDistributionFile::append(data,"CCMDF.out");
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
        x->factor[b] += step_size*(2*gsl_rng_uniform(r)-1);
    }
}

void P1(void *) {
    ORSA_DEBUG("print here...");
}



#endif // _VESTA_INTERIOR_ANALYTIC_H_
