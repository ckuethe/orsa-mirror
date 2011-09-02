#ifndef _VESTA_INTERIOR_ANALYTIC_H_
#define _VESTA_INTERIOR_ANALYTIC_H_

#include <orsa/util.h>
#include <orsa/statistic.h>

#include <orsaPDS/RadioScienceGravity.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_siman.h>

#include "CubicChebyshevMassDistribution.h"

#include <orsaUtil/adaptiveInterval.h>
#include <orsaUtil/adaptiveMonteCarlo.h>

// GSL Simulated Annealing

/* how many points do we try before stepping */      
#define N_TRIES 100 // 100 // 200             

/* how many iterations for each T? */
#define ITERS_FIXED_T 100 // 200 // 100 // 1000 // 

/* max step size in random walk */
#define STEP_SIZE 1.0           

/* Boltzmann constant */
#define K 1.0                   

/* initial temperature */
// #define T_INITIAL 0.008     
#define T_INITIAL 0.010         

/* damping factor for temperature */
#warning no damping??
#define MU_T 1.000 // 1.010 // 1.003      
#define T_MIN 1.0e-5 // 2.0e-6

gsl_siman_params_t params  = {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
                              K, T_INITIAL, MU_T, T_MIN};

class SIMAN_xp {
public:
    orsa::Cache<double> R0_plate;
    orsa::Cache<double> R0_gravity;
    orsa::Cache<double> bulkDensity;
    osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape;
    orsa::Cache<size_t> T_degree;
    orsa::Cache<size_t> T_size;
    gsl_vector * cT0;
    gsl_vector * * uK;
    orsa::Cache<size_t> uK_size;
    std::vector<double> factor;
};


void SIMAN_copy (void * source, void * dest) {
    SIMAN_xp * s = (SIMAN_xp *) source;
    SIMAN_xp * d = (SIMAN_xp *) dest;
    d->R0_plate            = s->R0_plate;
    d->R0_gravity          = s->R0_gravity;
    d->bulkDensity         = s->bulkDensity;
    d->randomPointsInShape = s->randomPointsInShape;
    d->T_degree            = s->T_degree;
    d->T_size              = s->T_size;
    d->cT0                 = s->cT0;
    d->uK                  = &(s->uK[0]);
    d->uK_size             = s->uK_size;
    d->factor              = s->factor;
}

void * SIMAN_copy_construct (void * xp) {
    SIMAN_xp * d = new SIMAN_xp;
    SIMAN_copy(xp,d);
    return d;
}

void SIMAN_destroy (void * xp) {
    delete (SIMAN_xp *) xp;
}

double E1(void * xp) {
    
    SIMAN_xp * x = (SIMAN_xp *) xp;
    
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
    
    osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
        new CubicChebyshevMassDistribution(coeff,x->bulkDensity,x->R0_plate);
    
    x->randomPointsInShape->updateMassDistribution(massDistribution);
    
    orsa::Vector v;
    double density;
    osg::ref_ptr< orsa::Statistic<double> > stat = new orsa::Statistic<double>;
    orsa::Cache<double> minDensity, maxDensity;
    x->randomPointsInShape->reset();
    while (x->randomPointsInShape->get(v,density)) { 
        stat->insert(density);
        minDensity.setIfSmaller(density);
        maxDensity.setIfLarger(density);
    }
    
#warning find a better way to compute average density? (based on simplexIntegral and not on randomPointsInShape)
    
    ORSA_DEBUG("[density] min: %+6.2f max: %+6.2f avg: %+6.2f [g/cm^3]",
               orsa::FromUnits(orsa::FromUnits(stat->min(),orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
               orsa::FromUnits(orsa::FromUnits(stat->max(),orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
               orsa::FromUnits(orsa::FromUnits(stat->average(),orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
    
    /* 
       {
       // quick output
       char filename[1024];
       sprintf(filename,"quickProfile_%+.6f_%d.dat",orsa::FromUnits(orsa::FromUnits(stat->min(),orsa::Unit::GRAM,-1),orsa::Unit::CM,3),(*orsa::GlobalRNG::randomSeed));
       ORSA_DEBUG("writing file [%s]",filename);
       FILE * fp = fopen(filename,"w");
       double PP = -285.0;
       while (PP < 285.0) {
       v = orsa::Vector(orsa::FromUnits(PP,orsa::Unit::KM),0,0);
       density = massDistribution->density(v);
       gmp_fprintf(fp,"%g %g\n",PP,orsa::FromUnits(orsa::FromUnits(density,orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
       PP += 1.0;
       }
       fclose(fp);
       }
    */

    if (x->bulkDensity*stat->min() >= 0.0) {
        // another quick output...
#warning pass filename as parameter...
        CubicChebyshevMassDistributionFile::CCMDF_data data;
        data.minDensity = stat->min();
        data.maxDensity = stat->max();
        data.deltaDensity = (stat->max()-stat->min());
        data.densityScale = x->bulkDensity;
        data.R0 = x->R0_plate;
        data.coeff = coeff;
        CubicChebyshevMassDistributionFile::append(data,"CCMDF.out");
    }
    
    // first approach: maximize the minimum density
    // return -minDensity;
    // alternative: minimize density range
    // return (stat->max()-stat->min());
    // more versions
    // return std::max(0.0,-minDensity);
    // return std::max(0.0,-minDensity)*(maxDensity-minDensity);
    // return -minDensity*(maxDensity-minDensity);
    return maxDensity-minDensity;
    // return -minDensity*(maxDensity-minDensity);
    // return minDensity*(maxDensity-minDensity);
    
    // maximum amplitude!
    // return minDensity-maxDensity;
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
        x->factor[b] += step_size*(2*gsl_rng_uniform(r)-1)/x->uK_size;
    }
}

void P1(void *) {
    ORSA_DEBUG("print here...");
}



#endif // _VESTA_INTERIOR_ANALYTIC_H_
