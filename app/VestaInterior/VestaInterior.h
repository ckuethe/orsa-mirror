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
    CubicChebyshevMassDistribution(const CoefficientType & coefficient) : orsa::MassDistribution(), coeff(coefficient) { }
protected:
    ~CubicChebyshevMassDistribution() { }
public:
    double density(const orsa::Vector & p) const {
        const size_t order = coeff.size()-1;
        double density = 0.0;
        std::vector<double> Tx, Ty, Tz;
        orsa::ChebyshevT(Tx,order,p.getX());
        orsa::ChebyshevT(Ty,order,p.getY());
        orsa::ChebyshevT(Tz,order,p.getZ());
        for (size_t i=0; i<=order; ++i) {
            for (size_t j=0; j<=order-i; ++j) {
                for (size_t k=0; k<=order-i-j; ++k) {
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
};

typedef Entry AdaptiveIntervalTemplateType;

typedef orsaUtil::AdaptiveIntervalElement<AdaptiveIntervalTemplateType> AdaptiveIntervalElementType;

class AdaptiveIntervalType : public orsaUtil::AdaptiveInterval<AdaptiveIntervalTemplateType> {
public:
    AdaptiveIntervalType(const double & min,
                         const double & max,
                         const double & confidenceLevel,
                         const double & initialThresholdLevel,
                         const double & targetThresholdLevel,
                         const size_t & targetSamples) :
        orsaUtil::AdaptiveInterval<AdaptiveIntervalTemplateType> (min,max,confidenceLevel,initialThresholdLevel,targetThresholdLevel,targetSamples)
        { }
public:
    void updateLevel(const AdaptiveIntervalElementType & e) const {
        
    }    
};






#endif // _VESTA_INTERIOR_H_
