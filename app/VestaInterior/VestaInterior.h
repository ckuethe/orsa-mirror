#ifndef _VESTA_INTERIOR_H_
#define _VESTA_INTERIOR_H_

#include <orsa/chebyshev.h>
#include <orsa/massDistribution.h>

#include <orsaUtil/adaptiveInterval.h>
#include <orsaUtil/adaptiveMonteCarlo.h>

#include <orsaPDS/RadioScienceGravity.h>

// Interior model part

class CubicChebyshevMassDistribution : public orsa::MassDistribution {
public:
    typedef std::vector< std::vector< std::vector<double> > > CoefficientType;
protected:
    CoefficientType coeff;
public:
    static void resize(CoefficientType & coeff, const size_t & degree) {
        coeff.resize(degree+1);
        for (size_t i=0; i<=degree; ++i) {
            coeff[i].resize(degree+1-i);
            for (size_t j=0; j<=degree-i; ++j) {
                coeff[i][j].resize(degree+1-i-j);
            }
        }
    }
public:
    CubicChebyshevMassDistribution(const CoefficientType & coefficient) : orsa::MassDistribution(), coeff(coefficient) { }
protected:
    ~CubicChebyshevMassDistribution() { }
public:
    double density(const orsa::Vector & p) const {
        const size_t degree = coeff.size()-1;
        double density = 0.0;
        std::vector<double> Tx, Ty, Tz;
        orsa::ChebyshevT(Tx,degree,p.getX());
        orsa::ChebyshevT(Ty,degree,p.getY());
        orsa::ChebyshevT(Tz,degree,p.getZ());
        for (size_t i=0; i<=degree; ++i) {
            for (size_t j=0; j<=degree-i; ++j) {
                for (size_t k=0; k<=degree-i-j; ++k) {
                    density += coeff[i][j][k]*Tx[i]*Ty[j]*Tz[k];
                }
            }
        }
        return density;
    }
};

// AdaptiveInterval part

class Entry : public osg::Referenced {
public:
    Entry() : osg::Referenced() { }
protected:
    ~Entry() { }
public:
    CubicChebyshevMassDistribution::CoefficientType coeff;
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
                         const size_t & targetSamples) :
        orsaUtil::AdaptiveInterval<AdaptiveIntervalTemplateType> (min,max,residualProbability,initialThresholdLevel,targetThresholdLevel,targetSamples)
        { }
public:
    void updateLevel(const AdaptiveIntervalElementType & e) const {
        if (e.level.isSet()) return;
        double chisq=0.0;
        
        // ...
        
        
        e.level = chisq;
    }    
};

typedef std::vector< osg::ref_ptr<AdaptiveIntervalType> > AdaptiveIntervalVector;

class AdaptiveMonteCarloType : public orsaUtil::AdaptiveMonteCarlo<AdaptiveIntervalVector> {
public:
    AdaptiveMonteCarloType() : orsaUtil::AdaptiveMonteCarlo<AdaptiveIntervalVector>() { }
protected:
    ~AdaptiveMonteCarloType() { }
    
    
};

#endif // _VESTA_INTERIOR_H_
