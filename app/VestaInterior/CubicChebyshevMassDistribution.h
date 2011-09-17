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
    const double densityScale;
public:
    static size_t totalSize(const size_t & degree);
public:
    static void resize(CoefficientType & coeff, const size_t & degree);
protected:
    typedef std::vector< std::vector< std::vector<size_t> > > IndexTableType;
    static IndexTableType indexTable;
public:
    static size_t index(const size_t & nx, const size_t & ny, const size_t & nz);
    static void triIndex(size_t & nx, size_t & ny, size_t & nz, const size_t & index);
protected:
    static void updateIndexTable(const size_t & requestedDegree);
public:
    CubicChebyshevMassDistribution(const CoefficientType & coefficient,
                                   const double & densityScale, // = bulk density if coefficients are relative
                                   const double & R0);
protected:
    virtual ~CubicChebyshevMassDistribution();
public:
    double density(const orsa::Vector & p) const;
};

// decompose a generic mass distribution into a cubic chebyshev mass distribuiton
// note: there are no tests on whether the points tested are inside or outside the body shape,
//       but that should not matter; what matters is that the new mass distribution
//       returns the same density as the input at any given point
CubicChebyshevMassDistribution * CubicChebyshevMassDistributionDecomposition(const orsa::MassDistribution * massDistribution,
                                                                             const size_t & degree,
                                                                             const double & densityScale,
                                                                             const double & R0);

class CubicChebyshevMassDistributionFile {
public:
    class CCMDF_data {
        // basic CCMD data + some auxil data useful to sort solutions
    public:
        double minDensity, maxDensity, deltaDensity;
        double penalty;
        double densityScale;
        double R0;
        size_t SH_degree;
        // size_T T_degree;
        CubicChebyshevMassDistribution::CoefficientType coeff;
    };
public:
    typedef CCMDF_data DataType;
    typedef std::deque<DataType> DataContainer;
public:
    static bool read(DataContainer & data, const std::string & fileName);
    static bool read(DataContainer & data, const std::string & fileName, const double & limitDeltaDensity);
public:
    static bool write(const DataContainer & data, const std::string & fileName);
    static bool write(const DataType & data, const std::string & fileName);
public:
    static bool append(const DataType & data, const std::string & fileName);
protected:
    static bool read(DataType & data, FILE * fp);
protected:
    static bool write(const DataType & data, FILE * fp);
};

#endif // CUBIC_CHEBYSHEV_MASS_DISTRIBUTION_H
