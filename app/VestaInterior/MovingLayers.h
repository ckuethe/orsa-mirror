#ifndef _MOVING_LAYERS_H_
#define _MOVING_LAYERS_H_

#include <orsa/util.h>
#include <orsa/statistic.h>

#include <orsaPDS/RadioScienceGravity.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_siman.h>

#include "CubicChebyshevMassDistribution.h"
#include "penalty.h"

#include "simplex.h"

typedef dd_real simplex_T;

// utils

// GSL Simulated Annealing

/* how many points do we try before stepping */      
#define N_TRIES 100 // 100 // 200             

/* how many iterations for each T? */
#define ITERS_FIXED_T 100 // 200 // 100 // 1000 // 

/* max step size in random walk */
#define STEP_SIZE 10000.0       

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
    std::vector<orsa::Vector> rv;
    orsa::Cache<size_t> SH_degree;
    orsa::Cache<size_t> T_degree;
    orsa::Cache<size_t> T_size;
    // gsl_vector * cT0;
    gsl_vector * * uK;
    orsa::Cache<size_t> uK_size;
    // std::vector<double> factor;
    orsa::Cache<double> minimumDensity;
    orsa::Cache<double> maximumDensity;
    orsa::Cache<double> penaltyThreshold;
    osg::ref_ptr<const LayerData> layerData;
    gsl_matrix * pseudoInvA;
    orsa::Vector sampled_CM;
    osg::ref_ptr< SimplexIntegration<simplex_T> > si;
    size_t M;
    size_t N;
    osg::ref_ptr<orsaPDS::RadioScienceGravityData> gravityData;
    gsl_vector * pds_coeff;
    std::vector< std::vector<mpf_class> > uniformShape_norm_C;
    std::vector< std::vector<mpf_class> > uniformShape_norm_S;
    // keep field entries in sync with copy function below!
};

void SIMAN_copy (void * source, void * dest) {
    SIMAN_xp * s = (SIMAN_xp *) source;
    SIMAN_xp * d = (SIMAN_xp *) dest;
    d->R0_plate            = s->R0_plate;
    d->R0_gravity          = s->R0_gravity;
    d->bulkDensity         = s->bulkDensity;
    d->rv                  = s->rv;
    d->SH_degree           = s->SH_degree;
    d->T_degree            = s->T_degree;
    d->T_size              = s->T_size;
    // d->cT0                 = s->cT0;
    d->uK                  = &(s->uK[0]);
    d->uK_size             = s->uK_size;
    // d->factor              = s->factor;
    d->minimumDensity      = s->minimumDensity;
    d->maximumDensity      = s->maximumDensity;
    d->penaltyThreshold    = s->penaltyThreshold;
    d->layerData           = s->layerData;
    d->pseudoInvA          = s->pseudoInvA;
    d->sampled_CM          = s->sampled_CM;
    d->si                  = s->si;
    d->M                   = s->M;
    d->N                   = s->N;
    d->gravityData         = s->gravityData;
    d->pds_coeff           = s->pds_coeff;
    d->uniformShape_norm_C = s->uniformShape_norm_C;
    d->uniformShape_norm_S = s->uniformShape_norm_S;
}

void * SIMAN_copy_construct (void * xp) {
    SIMAN_xp * d = new SIMAN_xp;
    SIMAN_copy(xp,d);
    return d;
}

void SIMAN_destroy (void * xp) {
    delete (SIMAN_xp *) xp;
}

double E1(void * xp);

double M1(void * xp, void * yp) {
    SIMAN_xp * x = (SIMAN_xp *) xp;
    SIMAN_xp * y = (SIMAN_xp *) yp;

#warning can improve this... use same scaling coefficients as in S1

    if (x->layerData->shLayerVector.size() != y->layerData->shLayerVector.size()) {
        ORSA_DEBUG("problems...");
    }
    
    double distance = 0.0;
    for (size_t k=0; k<x->layerData->shLayerVector.size(); ++k) {
        const osg::ref_ptr<LayerData::SHLayer> x_shL = x->layerData->shLayerVector[k];
        const osg::ref_ptr<LayerData::SHLayer> y_shL = y->layerData->shLayerVector[k];
        for (size_t l=0; l<x_shL->norm_A.size(); ++l) {
            for (size_t m=0; m<=l; ++m) {
                distance += orsa::square(x_shL->norm_A[l][m]-y_shL->norm_A[l][m]);
                if (m!=0) distance += orsa::square(x_shL->norm_B[l][m]-y_shL->norm_B[l][m]);
            }
        }
    }
    distance = sqrt(distance);
    return distance;
}

void S1(const gsl_rng * r, void * xp, double step_size) {
    
    ORSA_DEBUG("step size: %g",step_size);
    
    SIMAN_xp * x = (SIMAN_xp *) xp;
    
    LayerData::SHLayerVectorType newSHLayerVector;
    
    // modify norm_A and normB coefficients for each SH layer
    
    for (size_t k=0; k<x->layerData->shLayerVector.size(); ++k) {

        const osg::ref_ptr<LayerData::SHLayer> shL = x->layerData->shLayerVector[k];

        const double originalVolume = shL->volume();
        
        LayerData::SHLayer::SHcoeff new_norm_A = shL->norm_A;
        LayerData::SHLayer::SHcoeff new_norm_B = shL->norm_B;

        for (size_t l=0; l<new_norm_A.size(); ++l) {
            for (size_t m=0; m<=l; ++m) {
                new_norm_A[l][m] += step_size*(2*gsl_rng_uniform(r)-1);
                if (m!=0) new_norm_B[l][m] += step_size*(2*gsl_rng_uniform(r)-1);
            }
        }
        
        // make sure the mass fractions are unchanged, by scaling all norm_A and norm_B coefficients by a volume correction factor^(1/3)
        
        osg::ref_ptr<LayerData::SHLayer> tmpSHLayer = new LayerData::SHLayer(shL->excessDensity,
                                                                             new_norm_A,
                                                                             new_norm_B,
                                                                             shL->v0);
        const double tmpVolume = tmpSHLayer->volume();
        
        const double volumeCorrectionFactor = cbrt(originalVolume/tmpVolume);
        
        for (size_t l=0; l<new_norm_A.size(); ++l) {
            for (size_t m=0; m<=l; ++m) {
                new_norm_A[l][m] *= volumeCorrectionFactor;
                if (m!=0) new_norm_B[l][m] *= volumeCorrectionFactor;
            }
        }
        
        newSHLayerVector.push_back(new LayerData::SHLayer(shL->excessDensity,
                                                          new_norm_A,
                                                          new_norm_B,
                                                          shL->v0));
    }
    
    const LayerData::EllipsoidLayerVectorType oldEllipsoidLayerVector = x->layerData->ellipsoidLayerVector;
    x->layerData = new LayerData(oldEllipsoidLayerVector,
                                 newSHLayerVector);
    
}

void P1(void *) {
    ORSA_DEBUG("print here...");
}



#endif // _MOVING_LAYERS_H_
