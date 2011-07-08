#ifndef CUBIC_CHEBYSHEV_MASS_DISTRIBUTION_H
#define CUBIC_CHEBYSHEV_MASS_DISTRIBUTION_H

#include <orsa/chebyshev.h>
#include <orsa/massDistribution.h>

class CubicChebyshevMassDistribution : public orsa::MassDistribution {
public:
    typedef std::vector< std::vector< std::vector<double> > > CoefficientType;
protected:
    const CoefficientType coeff;
    const double oneOverR0;
public:
    static size_t totalSize(const size_t & degree);
public:
    static void resize(CoefficientType & coeff, const size_t & degree);
protected:
    static std::vector< std::vector< std::vector<size_t> > > indexTable;
public:
    static size_t index(const size_t & nx, const size_t & ny, const size_t & nz);
public:
    CubicChebyshevMassDistribution(const CoefficientType & coefficient,
                                   const double & R0);
protected:
    virtual ~CubicChebyshevMassDistribution();
public:
    double density(const orsa::Vector & p) const;
};

#endif // CUBIC_CHEBYSHEV_MASS_DISTRIBUTION_H
