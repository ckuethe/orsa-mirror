#ifndef CUBIC_CHEBYSHEV_MASS_DISTRIBUTION_H
#define CUBIC_CHEBYSHEV_MASS_DISTRIBUTION_H

#include <orsa/cache.h>
#include <orsa/chebyshev.h>
#include <orsa/massDistribution.h>
#include <orsa/legendre.h>
#include <orsa/util.h>

#include <string>

// note: most SH code is from SH2ijk.h, should merge into a single new header file

// "wedding cake" layers, with a base density, and excess densities on smaller and smaller volumes contaning each other
class LayerData : public osg::Referenced {
public:
    class EllipsoidLayer : public osg::Referenced {
    public:
        EllipsoidLayer(const double & excessDensity_,
                       const double & a_,
                       const double & b_,
                       const double & c_,
                       const orsa::Vector & v0_,
                       const orsa::Matrix & rot_) :
            osg::Referenced(true),
            excessDensity(excessDensity_),
            a(a_),
            b(b_),
            c(c_),
            v0(v0_),
            rot(rot_),
            inv_rot(orsa::Matrix::inverted(rot)),
            _am2(1.0/(a*a)),
            _bm2(1.0/(b*b)),
            _cm2(1.0/(c*c)),
            _volume(4.0/3.0*orsa::pi()*a*b*c),
            _excessMass(_volume*excessDensity) {
            // ORSA_DEBUG("excessDensity = %g",excessDensity);
            // ORSA_DEBUG("a,b,c = %g %g %g",a,b,c);
            // orsa::print(v0);
            // orsa::print(rot);
        }                     
    protected:
        virtual ~EllipsoidLayer() { }
    public: /* input */
        const double excessDensity;
        const double a,b,c; // semi-axes along x,y,z
        const orsa::Vector v0; // center of ellipsoid
        const orsa::Matrix rot; // rotation of ellipsoid
        const orsa::Matrix inv_rot; // inverse of rotation
    protected: /* derived */
        const double _am2,_bm2,_cm2;
        const double _volume;
        const double _excessMass;
    public:
        double am2() const { return _am2; }
        double bm2() const { return _bm2; }
        double cm2() const { return _cm2; }
        double volume() const { return _volume; }
        double excessMass() const { return _excessMass; }
    public:
        bool containsPoint(const orsa::Vector & p) const {
            // const orsa::Vector dp = p-v0;
            const orsa::Vector dp = inv_rot*(p-v0);
            return (orsa::square(dp.getX())*am2() +
                    orsa::square(dp.getY())*bm2() +
                    orsa::square(dp.getZ())*cm2() <= 1.0);
        } 
    public:
        /* bool containsLayer(const EllipsoidLayer * layer) const {
           #warning write better method...
           // simple one; a better one should consided difference in v0 and actual a,b,c values (tricky in particular cases...)
           return (volume() > layer->volume());
           }
        */
        // #warning errors should be generated in the code using the Layers if layer A is not inside layer B, and layer B is not inside layer A (i.e., they are crossing each other)
    };
public:
    typedef std::vector< osg::ref_ptr<EllipsoidLayer> > EllipsoidLayerVectorType;
    const EllipsoidLayerVectorType ellipsoidLayerVector;
    
public:
    // spherical harmonics layers
    class SHLayer : public osg::Referenced {
    public:  
        typedef std::vector< std::vector<double> > SHcoeff;
    public:
        SHLayer(const double & excessDensity_,
                const SHcoeff & norm_A_,
                const SHcoeff & norm_B_,
                const orsa::Vector & v0_,
                const orsa::Matrix & rot_) :
            osg::Referenced(true),
            excessDensity(excessDensity_),
            norm_A(norm_A_),
            norm_B(norm_B_),
            v0(v0_),
            rot(rot_),
            inv_rot(orsa::Matrix::inverted(rot))
            { }
    protected:
        virtual ~SHLayer() { }
    public: /* input */
        const double excessDensity;
        const SHcoeff norm_A, norm_B;
        const orsa::Vector v0; // center 
        const orsa::Matrix rot; // rotation
        const orsa::Matrix inv_rot; // inverse of rotation
    public:
        // derived
        mutable orsa::Cache<double> volume_, excessMass_;
        double volume() const;
        double excessMass() const;
    public:
        std::string MD5() const;
    protected:
        // norm_coeff = normalization_factor * coeff
        static double normalization_factor(const size_t & l,
                                           const size_t & m) {
            return orsa::normalization_sphericalHarmonicsToNormalizedSphericalHarmonics(l,m).get_d();
        }
    protected:
        // utility function, pure SH, v0 and rot NOT included
        double radius(const double & theta,
                      const double & phi) const {
            const double c_theta = cos(theta);
            double r=0;
            if (norm_A.size() != norm_B.size()) {
                ORSA_DEBUG("problems...");
            }
            for (size_t l=0; l<norm_A.size(); ++l) {
                if (norm_A[l].size() != norm_B[l].size()) {
                    ORSA_DEBUG("problems...");
                }
                for (size_t m=0; m<=l; ++m) {
                    if ( (norm_A[l][m] != 0.0) ||
                         (norm_B[l][m] != 0.0) ) {
                        const double Nf = normalization_factor(l,m);
                        const double Plm = orsa::LegendreP(l,m,c_theta).get_d()/Nf;
                        // ORSA_DEBUG("LegendreP[%i][%i](%g) = %g",l,m,c_theta,Plm);
                        r += Plm*(norm_A[l][m]*cos(m*phi));
                        if (m!=0) {
                            r += Plm*(norm_B[l][m]*sin(m*phi));
                        }   
                    }
                }
            }
            return r;
        }
        
    public:
        bool containsPoint(const orsa::Vector & p) const {
            // const orsa::Vector dp = p-v0;
            const orsa::Vector dp = inv_rot*(p-v0);
            const orsa::Vector  u = dp.normalized();
            const double theta = acos(u.getZ());
            const double phi   = atan2(u.getY(),u.getX());
            const double r     = radius(theta,phi);
            return (dp.lengthSquared() < r*r);
        }
    public:
        /* bool containsLayer(const EllipsoidLayer * layer) const {
           #warning write better method...
           } */
    };
public:
    typedef std::vector< osg::ref_ptr<SHLayer> > SHLayerVectorType;
    const SHLayerVectorType shLayerVector;
    
public:
    // check if layers are not crossing
    bool valid() {
        /* 
           const EllipsoidLayerVectorType & lv = ellipsoidLayerVector;
           for (unsigned int j=1; j<lv.size(); ++j) {
           for (unsigned int k=0; k<j; ++k) {
           if ( (!lv[k]->containsLayer(lv[j])) &&
           (!lv[j]->containsLayer(lv[k])) ) {
           return false;
           }
           }
           }
           #warning add check on SHLayer members...
        */
        return true;
    }
public:
    LayerData(const EllipsoidLayerVectorType & ellipsoidLayerVector_,
              const SHLayerVectorType & shLayerVector_) :
        osg::Referenced(true),
        ellipsoidLayerVector(ellipsoidLayerVector_),
        shLayerVector(shLayerVector_) { }
protected:
    virtual ~LayerData() { }
public:
    double density(const orsa::Vector & p) const {
        double density = 0.0;
        for (unsigned int k=0; k<ellipsoidLayerVector.size(); ++k) {
            if (ellipsoidLayerVector[k]->containsPoint(p)) {
                density += ellipsoidLayerVector[k]->excessDensity;
            }
        }
        for (unsigned int k=0; k<shLayerVector.size(); ++k) {
            if (shLayerVector[k]->containsPoint(p)) {
                density += shLayerVector[k]->excessDensity;
            }
        }
        return density;
    }
    double totalExcessMass() const {
        double M = 0.0;
        for (unsigned int k=0; k<ellipsoidLayerVector.size(); ++k) {
            M += ellipsoidLayerVector[k]->excessMass();
        }
        for (unsigned int k=0; k<shLayerVector.size(); ++k) {
            M += shLayerVector[k]->excessMass();
        }
        return M;
    }
};

// this class now also includes layerData data
class CubicChebyshevMassDistribution : public orsa::MassDistribution {
public:
    typedef std::vector< std::vector< std::vector<double> > > CoefficientType;
public:
    const CoefficientType coeff;
    // const double densityScale;
    const double oneOverR0;
    osg::ref_ptr<const LayerData> layerData;
public:
    static size_t totalSize(const size_t & degree);
public:
    static size_t degree(const CoefficientType & coeff);
    size_t degree() const;
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
                                   // const double & densityScale, // = bulk density if coefficients are relative
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
                                                                             // const double & densityScale,
                                                                             const double & R0,
                                                                             const LayerData * layerData = 0,
                                                                             const bool & decompose_layerData=false);

class CubicChebyshevMassDistributionFile {
public:
    class CCMDF_data {
        // basic CCMD data + some auxil data useful to sort solutions
    public:
        double minDensity, maxDensity, deltaDensity;
        double penalty;
        // double densityScale;
        double R0;
        size_t SH_degree;
        CubicChebyshevMassDistribution::CoefficientType coeff;
        osg::ref_ptr<const LayerData> layerData;
        std::string comment;
    public:
        void print() const;
    public:
        void clear() {
            minDensity = maxDensity = deltaDensity = 0.0;
            penalty = 0.0;
            // densityScale = 1.0;
            R0 = 1.0;
            SH_degree = 0;
            coeff.clear();
            layerData = 0;
            comment = "";
        }
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
                                              // data.densityScale,
                                              data.R0,
                                              data.layerData);
}

#endif // CUBIC_CHEBYSHEV_MASS_DISTRIBUTION_H
