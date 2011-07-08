#include "CubicChebyshevMassDistribution.h"

std::vector< std::vector< std::vector<size_t> > > CubicChebyshevMassDistribution::indexTable;

size_t CubicChebyshevMassDistribution::totalSize(const size_t & degree) {
    return (degree+1)*(degree+2)*(degree+3)/6;
}

void CubicChebyshevMassDistribution::resize(CoefficientType & coeff, const size_t & degree) {
    coeff.resize(degree+1);
    for (size_t i=0; i<=degree; ++i) {
        coeff[i].resize(degree+1-i);
        for (size_t j=0; j<=degree-i; ++j) {
            coeff[i][j].resize(degree+1-i-j);
        }
    }
}

size_t CubicChebyshevMassDistribution::index(const size_t & nx, const size_t & ny, const size_t & nz) {
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

CubicChebyshevMassDistribution::CubicChebyshevMassDistribution(const CoefficientType & coefficient,
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

CubicChebyshevMassDistribution::~CubicChebyshevMassDistribution() { }

double CubicChebyshevMassDistribution::density(const orsa::Vector & p) const {
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
