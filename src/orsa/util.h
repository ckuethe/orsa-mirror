#ifndef _ORSA_UTIL_
#define _ORSA_UTIL_

#include <string>

#include <orsa/body.h>
#include <orsa/bodygroup.h>
#include <orsa/cache.h>
#include <orsa/datetime.h>
#include <orsa/double.h>
#include <orsa/massDistribution.h>
#include <orsa/matrix.h>
#include <orsa/quaternion.h>
#include <orsa/shape.h>

#include <map>

#include <QHash>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace orsa {
  
    std::string & removeLeadingAndTrailingSpaces(std::string & s);
  
    std::string & removeAllSpaces(std::string & s);
  
    std::string & removeLeadingPlusSign(std::string & s);
  
    bool FrenetSerret(const orsa::Body * b,
                      orsa::BodyGroup  * bg,
                      const orsa::Time & t,
                      const orsa::Time & dt,
                      orsa::Vector & T,
                      orsa::Vector & N,
                      orsa::Vector & B);
  
    class IntPowCache {
    public:
        IntPowCache(const double & x) : _x(x) { }
    public:
        virtual ~IntPowCache() { }
    public:
        /* 
           const double get(const int p) const {
           if (!_data[p].isSet()) {
           // this could be more efficient...
           _data[p] = int_pow(_x,p);
           }
           return _data[p];
           }
        */
        //
        inline double get(const int p) const {      
            if (_data[p].isSet()) {
                return _data[p];
            } else {
                if (p > 0) {
                    _data[p] = get(p-1) * _x;
                } else if (p < 0) {
                    _data[p] = get(p+1) / _x;
                } else {
                    _data[p] = 1;
                }
                return _data[p];
            }
        }
    protected:
        const double _x;
    protected:
        // mutable QHash <int, orsa::Cache<double> > _data;
        mutable std::map <int, orsa::Cache<double> > _data;
    };
  
    bool eulerAnglesToMatrix(orsa::Matrix & m,
                             const double & psi,
                             const double & theta,
                             const double & phi);
  
    bool matrixToEulerAngles(double       & psi,
                             double       & theta,
                             double       & phi,
                             const orsa::Matrix & m);
  
    // wrapper
    /* 
       inline bool eulerAnglesToMatrix(orsa::Matrix                 & m,
       const RotationalBodyProperty * rotational) {
       if (rotational) {
       return eulerAnglesToMatrix(m, 
       rotational->getPsi(),
       rotational->getTheta(),
       rotational->getPhi());
       } else {
       // setting m should not be needed, but doesn't hurt
       m = orsa::Matrix::identity();
       return false;
       }
       }
    */
       
    // Find rotation matrix R that rotates unit vector a into unit vector b.
    // from: https://math.stackexchange.com/a/476311
    orsa::Matrix rotvec(const orsa::Vector & a, const orsa::Vector & b);
      
    orsa::Matrix QuaternionToMatrix (const orsa::Quaternion &);
  
    orsa::Quaternion MatrixToQuaternion (const orsa::Matrix &);
  
    orsa::Matrix localToGlobal(const orsa::Body       * b,
                               const orsa::BodyGroup  * bg,
                               const orsa::Time       & t);
  
    orsa::Matrix globalToLocal(const orsa::Body       * b,
                               const orsa::BodyGroup  * bg,
                               const orsa::Time       & t);
    
    // magnitude function
    // alpha = solar phase angle = angle Sun-Asteroid-Observer
    // G = slope parameter (G ~= 0.15)
    double P (const double & alpha, 
              const double & G);
    
    double apparentMagnitude(const double & H,
                             const double & G,
                             const double & phaseAngle,
                             const double & neo2obs,
                             const double & neo2sun);
    
    double absoluteMagnitude(const double & V,
                             const double & G,
                             const double & phaseAngle,
                             const double & neo2obs,
                             const double & neo2sun);
    
    // p = albedo, H = absolute magnitude
    double asteroidDiameter(const double & p, 
                            const double & H);
    
    void principalAxis(orsa::Matrix & genericToPrincipal,
                       orsa::Matrix & principalInertiaMatrix,
                       const orsa::Matrix & inertiaMatrix);
  
    // RNG class has two main purposes: 
    // - automatic memory handling; 
    // - state saving/restoring for checkpointing
    //
    // checkpoints: gsl_rng_fwrite and gsl_rng_fread (WARNING: binary files not portable across different architectures)
    // http://www.gnu.org/software/gsl/manual/html_node/Reading-and-writing-random-number-generator-state.html
    //
    // TODO: make this class thread safe?
    // 
    class RNG : public osg::Referenced  {
    public:
        RNG(int randomSeed_in) : 
            osg::Referenced(true),
            randomSeed(randomSeed_in) {
            commonInit();
        }
    protected:
        void commonInit() {
            rnd = ::gsl_rng_alloc(gsl_rng_gfsr4);
            ::gsl_rng_set(rnd,randomSeed);
        }
    protected:
        virtual ~RNG() {
            ::gsl_rng_free(rnd); 
        }
    public:
        int gsl_rng_fwrite(FILE * stream) const {
            return ::gsl_rng_fwrite(stream,rnd);
        }
    public:
        int gsl_rng_fread(FILE * stream) const {
            return ::gsl_rng_fread(stream,rnd);
        }
    public:
        double gsl_rng_uniform() const {
            return ::gsl_rng_uniform(rnd);
        }
    public:
        double gsl_rng_uniform_pos() const {
            return ::gsl_rng_uniform_pos(rnd);
        }
    public:
        unsigned long int gsl_rng_uniform_int(unsigned long int n) const {
            return ::gsl_rng_uniform_int(rnd,n);
        }
    public:
        void gsl_ran_dir_2d(double * x, double * y) const {
            ::gsl_ran_dir_2d(rnd,x,y);
        }
    public:
        void gsl_ran_dir_3d(double * x, double * y, double * z) const {
            ::gsl_ran_dir_3d(rnd,x,y,z);
        }
    public:
        double gsl_ran_binomial(double p, unsigned int n) const {
            return ::gsl_ran_binomial(rnd,p,n);
        }
    public:
        double gsl_ran_gaussian(double sigma) const {
            return ::gsl_ran_gaussian(rnd,sigma);
        }
    public:
        double gsl_ran_laplace(double a) const {
            return ::gsl_ran_laplace(rnd,a);
        }    
    public:
        unsigned int gsl_ran_poisson(double mu) const {
            return ::gsl_ran_poisson(rnd,mu);
        }
    public:
        void gsl_ran_dirichlet(size_t K, const double alpha[], double theta[]) const {
            return ::gsl_ran_dirichlet(rnd,K,alpha,theta);
        }
    public:
        const int randomSeed;
    protected:  
        gsl_rng * rnd;
    };
    
    /***/
    
    // A singletone class for the random number generator
    class GlobalRNG {   
    public:
        static GlobalRNG * instance() {
            if (_instance == 0) {
                _instance = new GlobalRNG;
            }
            return _instance;
        }
    protected:
        GlobalRNG() {
            if (!randomSeed.isSet()) {
                randomSeed=time(NULL)*getpid();
            }
            randomSeed.lock();
            ORSA_DEBUG("randomSeed: %i",(*randomSeed));
            rng_ = new orsa::RNG(randomSeed);
        }
    public:
        virtual ~GlobalRNG() {
            _instance = 0;
        }
    protected:
        static GlobalRNG * _instance;
    public:
        static orsa::Cache<int> randomSeed;
    protected:
        osg::ref_ptr<orsa::RNG> rng_;
    public:
        const orsa::RNG * rng() const { return rng_.get(); }
    public:
        void reset() const { randomSeed.unlock(); randomSeed.reset(); _instance=0; }
    };
    
    /***/
    
    class RandomPointsInShape : public osg::Referenced {
    public:
        RandomPointsInShape(const orsa::Shape * shape,
                            const orsa::MassDistribution * massDistribution,
                            const size_t & samplePoints,
                            const bool & localStoreVector=true);
    public:
        RandomPointsInShape(const orsa::TriShape * shape,
                            const orsa::MassDistribution * massDistribution,
                            const size_t & samplePoints,
                            const bool & localStoreVector=true);
    protected:
        virtual ~RandomPointsInShape() { } 
    public:
        osg::ref_ptr<const orsa::Shape> shape;
        osg::ref_ptr<const orsa::MassDistribution> md;
    public:
        const size_t size;
        const bool saveVector;
    public:
        size_t pointsInside() const { return numInside; }
    protected:
        size_t numInside;
    public:
        // if get(v) fails, you used all vectors
        // in order to use again all vectors, call reset()
        // bool get(orsa::Vector & v) const;
        bool get(orsa::Vector & v, double & density) const;
        bool get(orsa::Vector & v) const;
    public:
        void reset() const { 
            rng = new orsa::RNG(randomSeed);
            counter = 0; 
        }
    public:
        // keep list of vectors, just update the density
        void updateMassDistribution(const orsa::MassDistribution * massDistribution);
    protected:
        mutable size_t counter;
        mutable osg::ref_ptr<orsa::RNG> rng;
        static const int maxRandomSeed;
        const int randomSeed;
    protected:
        static orsa::Vector __randomVectorUtil(const orsa::RNG * rng,
                                               const Box & boundingBox);
    protected:
        std::vector<bool> in;
        std::vector< orsa::Cache<double> > density; 
        std::vector<orsa::Vector> vec;
    public:
        virtual RandomPointsInShape * clone() const {
            return new RandomPointsInShape(*this);
        }
    };
    
    //! this is the normalization coefficient needed to use the
    //! non-normalized spherical harmonics coefficients
    //! directly in the potential expansion formula
    inline mpf_class normalization_integralToSphericalHarmonics(const size_t & l, const size_t & m) {
        return mpf_class(mpq_class((2-orsa::kronecker(0,m))*orsa::factorial(l-m),
                                   (orsa::factorial(l+m))));
    }
    //! this converts SH coeff to normalized SH coeff
    inline mpf_class normalization_sphericalHarmonicsToNormalizedSphericalHarmonics(const size_t & l, const size_t & m) {
        return sqrt(mpf_class(mpq_class((orsa::factorial(l+m)),
                                        (2-orsa::kronecker(0,m))*(2*l+1)*orsa::factorial(l-m))));
    }
    //! this converts integrated coeff to normalized coeff, and is the product of the above two
    //! useful when converting from/to paul moments
    inline mpf_class normalization_integralToNormalizedSphericalHarmonics(const size_t & l, const size_t & m) {
        return sqrt(mpf_class(mpq_class((2-orsa::kronecker(m,0))*orsa::factorial(l-m), 
                                        (2*l+1)*orsa::factorial(l+m))));
    }
    
    double volume(const orsa::RandomPointsInShape * randomPointsInShape);
    
    double mass(const orsa::RandomPointsInShape * randomPointsInShape);
    
    orsa::Vector centerOfMass(const orsa::RandomPointsInShape * randomPointsInShape);
    // const orsa::MassDistribution * massDistribution);
    
    // InertiaMatrix about the centerOfMass, along the PrincipalAxis
    /* void principalAxisAndInertiaMatrix(orsa::Matrix & principalAxis,
       orsa::Matrix & inertiaMatrix,
       const orsa::Vector & centerOfMass,
       const orsa::RandomPointsInShape * randomPointsInShape,
       const orsa::MassDistribution * massDistribution);
    */
    //
    // InertiaMatrix about the centerOfMass
    // local here is the principal axis system
    void diagonalizedInertiaMatrix(orsa::Matrix & shapeToLocal,
                                   orsa::Matrix & localToShape,
                                   orsa::Matrix & inertiaMatrix,
                                   const orsa::Vector & centerOfMass,
                                   const orsa::RandomPointsInShape * randomPointsInShape);
    // const orsa::MassDistribution * massDistribution);
    
    orsa::PaulMoment * computePaulMoment(const unsigned int order,
                                         const orsa::Matrix & shapeToLocal,
                                         const orsa::Matrix & localToShape,
                                         const orsa::Vector & centerOfMass,
                                         const orsa::RandomPointsInShape * randomPointsInShape);
    // const orsa::MassDistribution * massDistribution);
    
    // utility, to perform all the computations above...
    void bodyInertialComputations(double & volume,
                                  orsa::Vector & centerOfMass,
                                  orsa::Matrix & shapeToLocal,
                                  orsa::Matrix & localToShape,
                                  orsa::Matrix & inertiaMatrix,
                                  orsa::PaulMoment * * paulMoment,
                                  const unsigned int order,
                                  const orsa::Shape * shape,
                                  const orsa::MassDistribution * massDistribution,
                                  const unsigned int N,
                                  const bool & localStoreVector=true);

    std::string randomString(const int len);
    
} // namespace orsa

#endif // _ORSA_UTIL_
