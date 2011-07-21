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

// AdaptiveInterval part

class AuxiliaryData : public osg::Referenced {
protected:
    ~AuxiliaryData() { }
public:
    orsa::Cache<double> R0;
    osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape;
    // orsa::Cache<double> bulkDensity;
    // orsa::Cache<double> bulkDensity_gcm3;
    orsa::Cache<size_t> chebyshevDegree;
    orsa::Cache<size_t> T_size;
    gsl_vector * cT0;
    gsl_vector * * uK;
    orsa::Cache<size_t> intervalVectorSize;
};

class Entry : public osg::Referenced {
public:
    Entry() : osg::Referenced() { }
protected:
    ~Entry() { }
public:
    CubicChebyshevMassDistribution::CoefficientType coeff;
public:
    void resize(const size_t & degree) {
        CubicChebyshevMassDistribution::resize(coeff,degree);
    }
};

typedef Entry AdaptiveIntervalTemplateType;

typedef orsaUtil::AdaptiveIntervalElement<AdaptiveIntervalTemplateType> AdaptiveIntervalElementType;

class AdaptiveIntervalType : public orsaUtil::AdaptiveInterval<AdaptiveIntervalTemplateType> {
public:
    AdaptiveIntervalType(const double & min,
                         const double & max,
                         const double & residualProbability,
                         const double & initialThresholdLevel,
                         const double & targetThresholdLevel,
                         const size_t & targetSamples,
                         const AuxiliaryData * auxiliaryData) :
        orsaUtil::AdaptiveInterval<AdaptiveIntervalTemplateType> (min,max,residualProbability,initialThresholdLevel,targetThresholdLevel,targetSamples),
        aux(auxiliaryData)
        { }
protected:
    osg::ref_ptr<const AuxiliaryData> aux;
protected:
    // mutable osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape;
public:
    void updateLevel(const AdaptiveIntervalElementType & e) const {

        if (e.level.isSet()) return;
        
        osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
            new CubicChebyshevMassDistribution(e.data->coeff,aux->R0);
        
        aux->randomPointsInShape->updateMassDistribution(massDistribution);
        
        orsa::Vector v;
        double density;
        osg::ref_ptr< orsa::Statistic<double> > stat = new orsa::Statistic<double>;
        aux->randomPointsInShape->reset();
        while (aux->randomPointsInShape->get(v,density)) { 
            stat->insert(density);
        }
        
        // e.level = std::max(-stat->min(),0.0);
        e.level = -stat->min();
        
        if (e.level <= 0.0) ORSA_DEBUG("level: %+8.3f",(*e.level));
    }
};

typedef std::vector< osg::ref_ptr<AdaptiveIntervalType> > AdaptiveIntervalVector;

class AdaptiveMonteCarloType : public orsaUtil::AdaptiveMonteCarlo<AdaptiveIntervalVector> {
public:
    AdaptiveMonteCarloType(const AuxiliaryData * auxiliaryData) :
        orsaUtil::AdaptiveMonteCarlo<AdaptiveIntervalVector>(),
        aux(auxiliaryData)
        { }
protected:
    ~AdaptiveMonteCarloType() { }
protected:
    osg::ref_ptr<const AuxiliaryData> aux;
public:
    bool sample(ElementTypeVector & ev) const {
        ++iterCount;
        // ORSA_DEBUG("iter: %i",iterCount);
        for (size_t k=0; k<aux->intervalVectorSize; ++k) {
            ev[k].position = intervalVector[k]->sample();
            ev[k].position.lock();
        }
        gsl_vector * cT = gsl_vector_calloc(aux->T_size);
        gsl_vector_memcpy(cT,aux->cT0);
        for (size_t k=0; k<aux->intervalVectorSize; ++k) {
            for (size_t l=0; l<aux->T_size; ++l) {
                gsl_vector_set(cT,l,gsl_vector_get(cT,l)+ev[k].position*gsl_vector_get(aux->uK[k],l));    
            }
        }
        // debug
        /* for (size_t l=0; l<aux->T_size; ++l) {
           size_t i,j,k;
           CubicChebyshevMassDistribution::triIndex(i,j,k,l);
           ORSA_DEBUG("cT[%i][%i][%i] = %+9.6f",i,j,k,gsl_vector_get(cT,l));    
           }
        */
        for (unsigned int p=0; p<aux->intervalVectorSize; ++p) {
            ev[p].data->resize(aux->chebyshevDegree);
            for (size_t l=0; l<aux->T_size; ++l) {
                size_t i,j,k;
                CubicChebyshevMassDistribution::triIndex(i,j,k,l);
                ev[p].data->coeff[i][j][k] = gsl_vector_get(cT,l);
            }
        }
        // compute level on one element, and copy it on all remaining elements
        intervalVector[0]->updateLevel(ev[0]);
        for (size_t p=1; p<aux->intervalVectorSize; ++p) {
            ev[p].level = ev[0].level;
        }
        for (size_t p=0; p<aux->intervalVectorSize; ++p) {
            ORSA_DEBUG("PosLev %3d %+8.3f %+8.3f",p,(*ev[p].position),(*ev[p].level));
        }
        gsl_vector_free(cT);
        return true;
    }
};

// GSL Simulated Annealing

/* how many points do we try before stepping */
#define N_TRIES 200             

/* how many iterations for each T? */
#define ITERS_FIXED_T 1000

/* max step size in random walk */
#define STEP_SIZE 1.0            

/* Boltzmann constant */
#define K 1.0                   

/* initial temperature */
#define T_INITIAL 0.008         

/* damping factor for temperature */
#define MU_T 1.003              
#define T_MIN 2.0e-6

gsl_siman_params_t params 
= {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
   K, T_INITIAL, MU_T, T_MIN};

class SIMAN_xp {
public:
    orsa::Cache<double> R0;
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
    d->R0                  = s->R0;
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
        for (unsigned int j=0; j<=x->T_degree; ++j) {
            for (unsigned int k=0; k<=x->T_degree; ++k) {
                if (i+j+k<=x->T_degree) {
                    const size_t index = CubicChebyshevMassDistribution::index(i,j,k);
                    coeff[i][j][k] = gsl_vector_get(cT,index);
                }
            }            
        }
    }
    
    
    osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
        new CubicChebyshevMassDistribution(coeff,x->R0);
    
    x->randomPointsInShape->updateMassDistribution(massDistribution);
    
    orsa::Vector v;
    double density;
    orsa::Cache<double> minDensity;
    x->randomPointsInShape->reset();
    while (x->randomPointsInShape->get(v,density)) { 
        minDensity.setIfSmaller(density);
    }
    
    return (-minDensity);
}

double M1(void * xp, void * yp) {
    SIMAN_xp * x = (SIMAN_xp *) xp;
    SIMAN_xp * y = (SIMAN_xp *) yp;

    double distance = 0.0;
    for (size_t b=0; b<x->uK_size; ++b) {
        distance += orsa::square(y->factor[b]-x->factor[b]);
    }
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
