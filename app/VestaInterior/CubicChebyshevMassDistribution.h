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
    typedef std::vector< std::vector< std::vector<size_t> > > IndexTableType;
    static IndexTableType indexTable;
public:
    static size_t index(const size_t & nx, const size_t & ny, const size_t & nz);
    static void triIndex(size_t & nx, size_t & ny, size_t & nz, const size_t & index);
protected:
    static void updateIndexTable(const size_t & requestedDegree);
public:
    CubicChebyshevMassDistribution(const CoefficientType & coefficient,
                                   const double & R0);
protected:
    virtual ~CubicChebyshevMassDistribution();
public:
    double density(const orsa::Vector & p) const;
};

class CubicChebyshevMassDistributionFile {
public:
    class CCMDF_data {
        // basic CCMD data + some auxil data useful to sort solutions
    public:
        double minDensity, maxDensity, deltaDensity;
        CubicChebyshevMassDistribution::CoefficientType coeff;
        // same R0 as plate model...
    };
public:
    typedef std::list<CCMDF_data> DataType;
public:
    static bool read(DataType & data, const std::string & fileName);
public:
    static bool write(const DataType & data, const std::string & fileName);
public:
    static bool append(const CCMDF_data & data, const std::string & fileName);
protected:
    static bool read(CCMDF_data & data, FILE * fp);
protected:
    static bool write(const CCMDF_data & data, FILE * fp);
};

#endif // CUBIC_CHEBYSHEV_MASS_DISTRIBUTION_H
