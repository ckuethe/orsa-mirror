#ifndef CUBIC_CHEBYSHEV_MASS_DISTRIBUTION_H
#define CUBIC_CHEBYSHEV_MASS_DISTRIBUTION_H

#include <orsa/chebyshev.h>
#include <orsa/massDistribution.h>

// "wedding cake" layers, with a base density, and excess densities on smaller and smaller volumes contaning each other
class LayerData : public osg::Referenced {
public:
    class EllipsoidLayer : public osg::Referenced {
    public:
        EllipsoidLayer(const double & excessDensity_,
                       const double & a_,
                       const double & b_,
                       const double & c_,
                       const orsa::Vector & v0_) :
            osg::Referenced(true),
            excessDensity(excessDensity_),
            a(a_),
            b(b_),
            c(c_),
            v0(v0_),
            am2(1.0/(a*a)),
            bm2(1.0/(b*b)),
            cm2(1.0/(c*c)),
            volume(4.0/3.0*orsa::pi()*a*b*c),
            excessMass(volume*excessDensity) { }
    protected:
        virtual ~EllipsoidLayer() { }
    public: /* input */
        const double excessDensity;
        const double a,b,c; // semi-axes along x,y,z
        const orsa::Vector v0; // center of ellipsoid
    public: /* derived */
        const double am2,bm2,cm2;
        const double volume;
        const double excessMass;
    public:
        bool containsPoint(const orsa::Vector & p) const {
            const orsa::Vector dp = p-v0;
            return (orsa::square(dp.getX())*am2 +
                    orsa::square(dp.getY())*bm2 +
                    orsa::square(dp.getZ())*cm2 <= 1.0);
        }
        bool containsLayer(const EllipsoidLayer * layer) const {
#warning write better method...
            // simple one; a better one should consided difference in v0 and actual a,b,c values (tricky in particular cases...)
            return (volume > layer->volume);
        }
#warning errors should be generated in the code using the Layers if layer A is not inside layer B, and layer B is not inside layer A (they are crossing each other)
    };
public:
    const double baseDensity;
public:
    typedef std::vector< osg::ref_ptr<EllipsoidLayer> > EllipsoidLayerVectorType;
    const EllipsoidLayerVectorType ellipsoidLayerVector;
public:
    // check if layers are not crossing
    bool valid() {
        const EllipsoidLayerVectorType & lv = ellipsoidLayerVector;
        for (unsigned int j=1; j<lv.size(); ++j) {
            for (unsigned int k=0; k<j; ++k) {
                if ( (!lv[k]->containsLayer(lv[j])) &&
                     (!lv[j]->containsLayer(lv[k])) ) {
                    return false;
                }
            }
        }
        return true;
    }
public:
    LayerData(const double & baseDensity_,
              const EllipsoidLayerVectorType & ellipsoidLayerVector_) :
        osg::Referenced(true),
        baseDensity(baseDensity_),
        ellipsoidLayerVector(ellipsoidLayerVector_) { }
protected:
    virtual ~LayerData() { }
public:
    // the density at a given point is baseDensity plus
    // the sum of all the excessDensities of all layers contining the point
    double density(const orsa::Vector & p) const {
        double density = baseDensity;
        const EllipsoidLayerVectorType & lv = ellipsoidLayerVector;
        for (unsigned int k=0; k<lv.size(); ++k) {
            if (lv[k]->containsPoint(p)) {
                density += lv[k]->excessDensity;
            }
        }
        return density;
    }
};

class CubicChebyshevMassDistribution : public orsa::MassDistribution {
public:
    typedef std::vector< std::vector< std::vector<double> > > CoefficientType;
public:
    const CoefficientType coeff;
    const double densityScale;
    const double oneOverR0;
    osg::ref_ptr<const LayerData> layerData;
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
                                   const double & R0,
                                   const LayerData * layerData = 0);
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
                                                                             const double & R0,
                                                                             const LayerData * layerData = 0);

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
        osg::ref_ptr<const LayerData> layerData;
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

inline CubicChebyshevMassDistribution * CCMD(const CubicChebyshevMassDistributionFile::CCMDF_data & data) {
    return new CubicChebyshevMassDistribution(data.coeff,
                                              data.densityScale,
                                              data.R0,
                                              data.layerData);
}

#endif // CUBIC_CHEBYSHEV_MASS_DISTRIBUTION_H
