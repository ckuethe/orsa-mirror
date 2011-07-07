#ifndef _VESTA_INTERIOR_H_
#define _VESTA_INTERIOR_H_

#include <orsa/chebyshev.h>
#include <orsa/massDistribution.h>
#include <orsa/util.h>

#include <orsaUtil/adaptiveInterval.h>
#include <orsaUtil/adaptiveMonteCarlo.h>

#include <orsaPDS/RadioScienceGravity.h>

#include <gsl/gsl_blas.h>

// Interior model part

class CubicChebyshevMassDistribution : public orsa::MassDistribution {
public:
    typedef std::vector< std::vector< std::vector<double> > > CoefficientType;
protected:
    const CoefficientType coeff;
    const double oneOverR0;
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
protected:
    static std::vector< std::vector< std::vector<size_t> > > indexTable;
public:
    static size_t index(const size_t & nx, const size_t & ny, const size_t & nz) {
        const size_t requestedDegree=nx+ny+nz;
        if (indexTable.size() < (requestedDegree+1)) {
            indexTable.resize(requestedDegree+1);
            for (size_t i=0; i<=requestedDegree; ++i) {
                indexTable[i].resize(requestedDegree+1-i);
                for (size_t j=0; j<=requestedDegree-i; ++j) {
                    indexTable[i][j].resize(requestedDegree+1-i-j);
                }
            }
            size_t idx=0;
            size_t degree=0;
            while (degree <= requestedDegree) {
                for (unsigned int i=0; i<=degree; ++i) {
                    for (unsigned int j=0; j<=degree; ++j) {
                        for (unsigned int k=0; k<=degree; ++k) {
                            if (i+j+k==degree) {
                                // ORSA_DEBUG("inserting %i-%i-%i  index: %i",i,j,k,idx);
                                indexTable[i][j][k] = idx++;
                            }
                        }
                    }
                }
                ++degree;
            }
        }
        return indexTable[nx][ny][nz];
    }
public:
    CubicChebyshevMassDistribution(const CoefficientType & coefficient,
                                   const double & R0) :
        orsa::MassDistribution(),
        coeff(coefficient),
        oneOverR0(1.0/R0) {
        /* const size_t degree = coeff.size()-1;
           for (unsigned int printDegree=0; printDegree<=degree; ++printDegree) {
           for (unsigned int i=0; i<=degree; ++i) {
           for (unsigned int j=0; j<=degree; ++j) {
           for (unsigned int k=0; k<=degree; ++k) {
           if (i+j+k==printDegree) {
           ORSA_DEBUG("coeff[%i][%i][%i] = %g",i,j,k,coeff[i][j][k]);
           }
           }
           }
           }
           }
        */
    }
protected:
    ~CubicChebyshevMassDistribution() { }
public:
    double density(const orsa::Vector & p) const {
        if (0) {
            // debug
            static size_t calls=0;
            ++calls;
            if (calls%1000==0) ORSA_DEBUG("calls: %i",calls);
        }
        if (coeff.size() == 0) return 0.0;
        const size_t degree = coeff.size()-1;
        std::vector<double> Tx, Ty, Tz;
        orsa::ChebyshevT(Tx,degree,p.getX()*oneOverR0);
        orsa::ChebyshevT(Ty,degree,p.getY()*oneOverR0);
        orsa::ChebyshevT(Tz,degree,p.getZ()*oneOverR0);
        double density = 0.0;
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

class AuxiliaryData : public osg::Referenced {
protected:
    ~AuxiliaryData() { }
public:
    osg::ref_ptr<orsa::Shape> shape;
    orsa::Cache<size_t> sphericalHarmonicDegree; // degree of the spherical harmonics expansion
    orsa::Cache<size_t> chebyshevDegree;
    orsa::Cache<double> R0;
    orsa::Cache<size_t> numSamplePoints;
    orsa::Cache<bool>   storeSamplePoints;
    orsa::Cache<size_t> intervalVectorSize;
    osg::ref_ptr<orsaPDS::RadioScienceGravityData> pds_data;
    orsa::Cache<size_t> pds_numberOfCoefficients;
    gsl_vector * pds_coeff;
    gsl_matrix * pds_inv_covm;
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
    mutable osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape;
public:
    void updateLevel(const AdaptiveIntervalElementType & e) const {
        if (e.level.isSet()) return;
        
        osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
            new CubicChebyshevMassDistribution(e.data->coeff,aux->R0);

        if (randomPointsInShape.get() == 0) {
            // first call only
            randomPointsInShape = new orsa::RandomPointsInShape(aux->shape,massDistribution,aux->numSamplePoints,aux->storeSamplePoints);
        } else {
            randomPointsInShape->updateMassDistribution(massDistribution);
        }
        
        // test: any point with negative mass?
        {
            orsa::Vector v;
            double density;
            randomPointsInShape->reset();
            while (randomPointsInShape->get(v,density)) { 
                if (density < 0.0) {
                    // ORSA_DEBUG("negative density...");
                    e.level = 1.e100;
                    return;
                }
            }
        }

        {
            // print out coeff
            const size_t degree = e.data->coeff.size()-1;
            for (unsigned int printDegree=0; printDegree<=degree; ++printDegree) {
                for (unsigned int i=0; i<=degree; ++i) {
                    for (unsigned int j=0; j<=degree; ++j) {
                        for (unsigned int k=0; k<=degree; ++k) {
                            if (i+j+k==printDegree) {
                                ORSA_DEBUG("coeff[%i][%i][%i] = %g",i,j,k,e.data->coeff[i][j][k]);
                            }
                        }
                    }
                }
            }
        }
        
        // const double volume = orsa::volume(randomPointsInShape);
        
        // check: all points inside shape have positive density?
        
        const double mass = orsa::mass(randomPointsInShape);

        const double volume = orsa::volume(randomPointsInShape);
        
        ORSA_DEBUG("bulk density: %g [g/cm^3]   volume: %g [km^3]",
                   orsa::FromUnits(orsa::FromUnits(mass/volume,orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
                   orsa::FromUnits(volume,orsa::Unit::KM,-3));
        
        const double GM = orsa::Unit::G()*mass;
        
        /* ORSA_DEBUG("mass: %g [kg] = %g [MSun]",
           orsa::FromUnits(mass,orsa::Unit::KG,-1),
           mass/orsaSolarSystem::Data::MSun());
           ORSA_DEBUG("GM: %g [km^3/s^2]",
           orsa::FromUnits(orsa::FromUnits(mass*orsa::Unit::G(),orsa::Unit::KM,-3),orsa::Unit::SECOND,2));
        */
        
        const orsa::Vector centerOfMass = orsa::centerOfMass(randomPointsInShape);
        // massDistribution);
        
        /* orsa::diagonalizedInertiaMatrix(shapeToLocal,
           localToShape,
           inertiaMatrix,
           centerOfMass,
           randomPointsInShape,
           massDistribution);
        */
        
        osg::ref_ptr<orsa::PaulMoment> paulMoment =
            orsa::computePaulMoment(aux->sphericalHarmonicDegree,
                                    orsa::Matrix::identity(), // shapeToLocal,
                                    orsa::Matrix::identity(), // localToShape,
                                    centerOfMass,
                                    randomPointsInShape);
        // massDistribution);
        
        std::vector< std::vector<double> > C, S, norm_C, norm_S;
        std::vector<double> J;
        orsa::convert(C, S, norm_C, norm_S, J,
                      paulMoment, 
                      aux->R0);
        
        /* for (size_t k=0; k<aux->pds_numberOfCoefficients; ++k) {
           ORSA_DEBUG("pds_coeff[%03i] = %g",k,gsl_vector_get(aux->pds_coeff,k));
           }
        */
        
        gsl_vector * vec_coeff = gsl_vector_alloc(aux->pds_numberOfCoefficients);
        orsa::Cache<size_t> minIndex, maxIndex;
        {
            size_t index;
            // first: GM
            index = aux->pds_data->index("GM");
            minIndex.setIfSmaller(index);
            maxIndex.setIfLarger(index);
            // ORSA_DEBUG("index: %i",index);
            ORSA_DEBUG("GM: %g",GM);
            ORSA_DEBUG("pds_coeff[%03i] = %g",index,gsl_vector_get(aux->pds_coeff,index));
            gsl_vector_set(vec_coeff,
                           index,
                           GM-gsl_vector_get(aux->pds_coeff,index));
            for (size_t l=2; l<=aux->sphericalHarmonicDegree; ++l) {
                for (size_t m=0; m<=l; ++m) {
                    index = aux->pds_data->index(orsaPDS::RadioScienceGravityData::keyC(l,m));
                    minIndex.setIfSmaller(index);
                    maxIndex.setIfLarger(index);
                    // ORSA_DEBUG("index: %i",index);
                    ORSA_DEBUG("norm_C[%03i][%03i] = %g",l,m,norm_C[l][m]);
                    ORSA_DEBUG("pds_coeff[%03i] = %g",index,gsl_vector_get(aux->pds_coeff,index));
                    gsl_vector_set(vec_coeff,
                                   index,
                                   norm_C[l][m]-gsl_vector_get(aux->pds_coeff,index));
                    if (m!=0) {
                        index = aux->pds_data->index(orsaPDS::RadioScienceGravityData::keyS(l,m));
                        minIndex.setIfSmaller(index);
                        maxIndex.setIfLarger(index);
                        // ORSA_DEBUG("index: %i",index);
                        ORSA_DEBUG("norm_S[%03i][%03i] = %g",l,m,norm_S[l][m]);
                        ORSA_DEBUG("pds_coeff[%03i] = %g",index,gsl_vector_get(aux->pds_coeff,index));
                        gsl_vector_set(vec_coeff,
                                       index,
                                       norm_S[l][m]-gsl_vector_get(aux->pds_coeff,index));
                    }
                }
            }
        }
        
        // ORSA_DEBUG("index range: %i -> %i [inclusive]",(*minIndex),(*maxIndex));
        
        if (minIndex != 0) {
            ORSA_DEBUG("this is unexpected...");
        }
        
        const size_t view_size = maxIndex+1-minIndex;
        gsl_vector_const_view vec_coeff_view = gsl_vector_const_subvector(vec_coeff,0,view_size);
        gsl_matrix_const_view inv_covm_view = gsl_matrix_const_submatrix(aux->pds_inv_covm,0,0,view_size,view_size);
        
        for (size_t k=0; k<view_size; ++k) {
            ORSA_DEBUG("vec_coeff[%03i] = %g",k,gsl_vector_get(vec_coeff,k));
        }
        
        
        gsl_vector * vec_tmp = gsl_vector_alloc(view_size);
        double chisq; 
        // gsl_blas_dgemv(CblasNoTrans,1.0,aux->pds_inv_covm,vec_coeff,0.0,vec_tmp);
        gsl_blas_dgemv(CblasNoTrans,1.0,&inv_covm_view.matrix,&vec_coeff_view.vector,0.0,vec_tmp);
        // gsl_blas_ddot (vec_coeff,vec_tmp,&chisq);
        gsl_blas_ddot(&vec_coeff_view.vector,vec_tmp,&chisq);
        
        e.level = chisq;
        
        ORSA_DEBUG("chisq = %g",chisq);
        
        gsl_vector_free(vec_coeff);
        gsl_vector_free(vec_tmp);
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
        // debug only
        // if (iterCount==3) exit(0);
        // #warning ^^^^^^^^^^^^^^^comment out line above!
        for (size_t k=0; k<aux->intervalVectorSize; ++k) {
            ev[k].position = intervalVector[k]->sample();
            ev[k].position.lock();
        }
        for (unsigned int p=0; p<aux->intervalVectorSize; ++p) {
            ev[p].data->resize(aux->chebyshevDegree);
            for (unsigned int i=0; i<=aux->chebyshevDegree; ++i) {
                for (unsigned int j=0; j<=aux->chebyshevDegree; ++j) {
                    for (unsigned int k=0; k<=aux->chebyshevDegree; ++k) {
                        if (i+j+k<=aux->chebyshevDegree) {
                            ev[p].data->coeff[i][j][k] = ev[CubicChebyshevMassDistribution::index(i,j,k)].position;
                            // ORSA_DEBUG("i: %i  j: %i  k: %i  index: %i",i,j,k,CubicChebyshevMassDistribution::index(i,j,k));
                        }
                    }
                }
            }
        }
        // compute level on one element, and copy it on all remaining elements
        intervalVector[0]->updateLevel(ev[0]);
        for (size_t p=0; p<aux->intervalVectorSize; ++p) {
            ev[p].level = ev[0].level;
        }
        return true;
    }
};

#endif // _VESTA_INTERIOR_H_
