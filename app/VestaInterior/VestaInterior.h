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
    static size_t totalSize(const size_t & degree) {
        return (degree+1)*(degree+2)*(degree+3)/6;
    }
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
    static size_t index(const size_t & nx, const size_t & ny, const size_t & nz) {
        const size_t requestedDegree=nx+ny+nz;
        if (indexTable.size() >= (requestedDegree+1)) {
            // ORSA_DEBUG("found: %i",indexTable[nx][ny][nz]);
            return indexTable[nx][ny][nz];
        } else {
            indexTable.resize(requestedDegree+1);
            for (size_t i=0; i<=requestedDegree; ++i) {
                indexTable[i].resize(requestedDegree+1-i);
                for (size_t j=0; j<=requestedDegree-i; ++j) {
                    indexTable[i][j].resize(requestedDegree+1-i-j);
                }
            }
            size_t idx=0;
            size_t degree=0;
            bool done=false;
            while (!done) {
                for (unsigned int i=0; i<=degree; ++i) {
                    for (unsigned int j=0; j<=degree; ++j) {
                        for (unsigned int k=0; k<=degree; ++k) {
                            if (i+j+k==degree) {
                                // ORSA_DEBUG("inserting %i-%i-%i  index: %i",i,j,k,idx);
                                indexTable[i][j][k] = idx;
                                if ((i==nx) && (j==ny) && (k==nz)) {
                                    done=true;
                                } else {
                                    ++idx;
                                }
                            }
                            if (done) break;
                        }
                        if (done) break;
                    }
                    if (done) break;
                }
                ++degree;
            }
            return idx;
        }
    }
protected:
    static std::vector< std::vector< std::vector<size_t> > > indexTable;
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
        
        osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
            new CubicChebyshevMassDistribution(e.data->coeff);
        
        
        /*

        orsa::Vector centerOfMass;
        orsa::Matrix shapeToLocal;
        orsa::Matrix localToShape;
        orsa::Matrix inertiaMatrix;
        osg::ref_ptr<orsa::PaulMoment> paulMoment;
        
        const unsigned int order = 4;
        const unsigned int N = 100000;
    
    // osg::ref_ptr<RandomPointsInShape> randomPointsInShape = new RandomPointsInShape(shape,N,randomSeed);
    osg::ref_ptr<RandomPointsInShape> randomPointsInShape = new RandomPointsInShape(shape,N);
    
    const double volume = orsa::volume(randomPointsInShape);
    
    centerOfMass = orsa::centerOfMass(randomPointsInShape,
                                      massDistribution);
    
    orsa::diagonalizedInertiaMatrix(shapeToLocal,
                                    localToShape,
                                    inertiaMatrix,
                                    centerOfMass,
                                    randomPointsInShape,
                                    massDistribution);
    
    paulMoment = orsa::computePaulMoment(order,
                                         shapeToLocal,
                                         localToShape,
                                         centerOfMass,
                                         randomPointsInShape,
                                         massDistribution);
    
	std::vector< std::vector<double> > C, S, norm_C, norm_S;
	std::vector<double> J;
	orsa::convert(C, S, norm_C, norm_S, J,
                  paulMoment, 
                  FromUnits(300,orsa::Unit::KM));
    
        */
        double chisq=0.0;
        
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
