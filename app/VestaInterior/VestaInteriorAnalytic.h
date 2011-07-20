#ifndef _VESTA_INTERIOR_ANALYTIC_H_
#define _VESTA_INTERIOR_ANALYTIC_H_

#include <orsa/util.h>
#include <orsa/statistic.h>

#include <orsaPDS/RadioScienceGravity.h>

#include <gsl/gsl_blas.h>

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
    orsa::Cache<double> bulkDensity;
    orsa::Cache<double> bulkDensity_gcm3;
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
            for (unsigned int i=0; i<=aux->chebyshevDegree; ++i) {
                for (unsigned int j=0; j<=aux->chebyshevDegree; ++j) {
                    for (unsigned int k=0; k<=aux->chebyshevDegree; ++k) {
                        if (i+j+k<=aux->chebyshevDegree) {
                            ev[p].data->coeff[i][j][k] = gsl_vector_get(cT,CubicChebyshevMassDistribution::index(i,j,k));
                        }
                    }
                }
            }
        }
        // compute level on one element, and copy it on all remaining elements
        intervalVector[0]->updateLevel(ev[0]);
        for (size_t p=1; p<aux->intervalVectorSize; ++p) {
            ev[p].level = ev[0].level;
        }
        gsl_vector_free(cT);
        return true;
    }
};

#endif // _VESTA_INTERIOR_ANALYTIC_H_
