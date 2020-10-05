#ifndef _ORSA_SHAPE_
#define _ORSA_SHAPE_

#include <osg/Referenced>

#include <orsa/box.h>
#include <orsa/cache.h>
#include <orsa/debug.h>
#include <orsa/double.h>
#include <orsa/vector.h>

#include <gsl/gsl_roots.h>

#include <vector>

namespace orsa {
  
    class Shape : public osg::Referenced {
    
    public:
        Shape() : Referenced(true) { }
    protected:
        virtual ~Shape() { }
    
    public:
        enum ShapeType {
            SHAPE_NONE,
            SHAPE_TRI,
            SHAPE_LATLON,
            SHAPE_ELLIPSOID
        };
    public:
        virtual ShapeType getType() const = 0;
    
        /* 
           public:
           enum NormalType {
           NORMAL_RADIAL,       // parallel to the vertex
           NORMAL_FACE,         // normal to the face
           NORMAL_FACE_AVERAGE, // averaged on all the faces containing a vertex
           NORMAL_DEFAULT=NORMAL_FACE_AVERAGE
           };
           public:
           const NormalType getNormalType() const {
           return _normal_type;
           }
           bool setNormalType(const NormalType type) {
           _normal_type = type;
           // _updateNormal();
           return true;  
           }
           protected:
           NormalType _normal_type;
           protected:
           // virtual bool _updateNormal() const = 0;
       
        */
    
        /* 
           public:
           virtual bool read() = 0;
        */
    
    public:
        virtual bool isInside(const Vector &) const = 0;
    
    public:
        virtual const Vector closestVertex(const Vector &) const = 0;
        
    public:
        virtual bool rayIntersection(orsa::Vector & intersectionPoint,
                                     orsa::Vector & normal,
                                     const orsa::Vector & P,
                                     const orsa::Vector & u,
                                     const bool fullLine = false) const = 0;
    public:
        virtual double volume() const = 0;
        
    public:
        const double & boundingRadius() const { 
            if (!_r_max.isSet()) {
                if (!_updateCache()) {
                    ORSA_DEBUG("problems encountered while updating cached values");
                }
            }
            return _r_max;
        }
    public:
        const orsa::Box & boundingBox() const { 
            if (!_boundingBox.isSet()) {
                if (!_updateCache()) {
                    ORSA_DEBUG("problems encountered while updating cached values");
                }
            }
            return _boundingBox;
        }
    public:
        const orsa::Box & symmetricBoundingBox() const { 
            if (!_symmetricBoundingBox.isSet()) {
                if (!_updateCache()) {
                    ORSA_DEBUG("problems encountered while updating cached values");
                }
            }
            return _symmetricBoundingBox;
        }
    protected:
        virtual bool _updateCache() const = 0;
    protected:
        mutable orsa::Cache<double> _r_min, _r_max;
    protected:    
        mutable orsa::Box _boundingBox;
    protected:    
        mutable orsa::Box _symmetricBoundingBox;
    };
    
    class TriShape : public orsa::Shape {
    public:
        class TriIndex {
        public:
            TriIndex(const unsigned int i, 
                     const unsigned int j, 
                     const unsigned int k) : 
                _i(i), 
                _j(j),
                _k(k) { }
        public:
            unsigned int i() const { return _i; }
        public:
            unsigned int j() const { return _j; }
        public:
            unsigned int k() const { return _k; }
        protected:
            unsigned int _i, _j, _k;
        };
    public:
        typedef std::vector<orsa::Vector> VertexVector;
        typedef std::vector<TriIndex> FaceVector;
    public:
        TriShape() : Shape() {
            _init();
        }
    public:
        TriShape(const VertexVector & v,
                 const FaceVector   & f) : 
            Shape(),
            _vertex(v),
            _face(f) {
            _init();
        }
    private:
        void _init() {
            _old_closest_vertex_index = 0; 
        }
    protected:
        ~TriShape() { }
    
    public:
        ShapeType getType() const {
            return SHAPE_TRI;
        }
    
    protected:
        bool _updateCache() const;
    
    public:
        bool isInside(const orsa::Vector &) const;
    protected:
        //! This method uses internal reference points.
        bool _isInside_useLineMethod(const Vector &) const;
    protected:
        //! This method assumes that the normals always point outward.
        bool _isInside_useNormalMethod(const Vector &) const;
    protected:    
        //
        bool _isInside_useFaceMethod(const Vector &) const;
    protected:
        //
        bool _isInside_usePointInTetrahedronMethod(const Vector &) const;
        
    public:
        const Vector closestVertex(const Vector &) const;   
        unsigned int closestVertexIndex(const Vector &) const;   
    protected:
        mutable unsigned int _old_closest_vertex_index;
    
    public:
        bool rayIntersection(orsa::Vector & intersectionPoint,
                             orsa::Vector & normal,
                             const orsa::Vector & P,
                             const orsa::Vector & u,
                             const bool fullLine = false) const;
    public:
        bool rayIntersection(orsa::Vector & intersectionPoint,
                             orsa::Vector & normal,
                             unsigned int & faceIndex,
                             const orsa::Vector & P,
                             const orsa::Vector & u,
                             const bool fullLine = false) const;
    public:
        double volume() const;
    protected:
        mutable orsa::Cache<double> _volume;
        
    protected:
        mutable orsa::Cache<double> _delta_min, _delta_max;
    
    protected:
        // reference points inside model, used by _isInside_useLineMethod
        mutable std::vector<Vector> _ref;
    
    protected:   
        VertexVector _vertex;
    public:
        const VertexVector & getVertexVector() const {
            return _vertex;
        }
    
    protected:    
        FaceVector _face;
    public:
        const FaceVector & getFaceVector() const {
            return _face;
        }  
        
    protected:
        mutable std::vector< std::list<unsigned int> > vertexInFace;
        
    public:
        //! This normal is the average of the normals of all the faces containing this vertex
        const Vector & _getVertexNormal(const unsigned int vertex_index) const;
    protected:    
        mutable std::vector<Vector> _vertex_normal; // size: [vertex]
    
    public:
        const Vector & _getFaceNormal(const unsigned int face_index) const;
    protected:  
        mutable std::vector<Vector> _face_normal;   // size: [face]
    
    public:
        double _getFaceArea(const unsigned int face_index) const;
    protected:  
        mutable std::vector<double> _face_area;   // size: [face]
        
    public:
        double surfaceArea() const;
    protected:
        mutable orsa::Cache<double> _area;
    
    public:
        typedef std::vector<double> AngleVector;
    public:
        //! lightSource and observerPosition are relative to the Shape
        //! delta is the angle, at the observer, between the center of the shape and the vertex
        bool vertexIlluminationAngles(const orsa::Vector & lightSource,
                                      const orsa::Vector & observerPosition,
                                      double             & phase,
                                      AngleVector        & i, 
                                      AngleVector        & e,
                                      AngleVector        & delta,
                                      const double       & deltaMax,
                                      const bool includeShadows = false) const;
    public:
        /* // this needs to be perfectioned before use
           bool faceIlluminationAngles(const orsa::Vector & lightSource,
           const orsa::Vector & observerPosition,
           double       & phase,
           AngleVector        & i, 
           AngleVector        & e,
           AngleVector        & delta,
           const double & deltaMax,
           const bool includeShadows = false) const;
        */
        
        // geodesic grid
        // generates the VertexVector and FaceVector to input in a new TriShape object
        // Nsub is the number of subdivisions
    public:
        static void GeodesicGrid(VertexVector & v,
                                 FaceVector   & f,
                                 const size_t & Nsub,
                                 const bool   & verbose=false);
        
    };
    
    class LatLonShape : public orsa::Shape {
    public:
        LatLonShape() : Shape() { }
    protected:
        ~LatLonShape() { }
    
    public:
        ShapeType getType() const {
            return SHAPE_LATLON;
        }
    
    protected:
        bool _updateCache() const;
    
    public:
        bool isInside(const orsa::Vector &) const;
        
    public:
        double volume() const { ORSA_DEBUG("code needed!"); return 0.0; }
        
    protected:
        // lat: -90 to +90, lon: 0 to 360
        double _lat_step_DEG;
        double _lon_step_DEG;
    public:
        const double & getLatStepDEG() const { return _lat_step_DEG; }
        const double & getLonStepDEG() const { return _lon_step_DEG; }
    
    public:
        typedef std::vector < std::vector< double > > RadiusTable;
    protected:
        RadiusTable _rt;
    public:
        const RadiusTable & getRadiusTable() const {
            return _rt;
        }
    };
  
    class EllipsoidShape : public orsa::Shape {
    public:
        //! a,b,c are the three semi-axis respectively relative to x,y,z  
        EllipsoidShape(const double & a,
                       const double & b,
                       const double & c) :
            Shape(),
            _a(a),
            _a2(a*a),
            _am2(1/_a2),
            _b(b),
            _b2(b*b),
            _bm2(1/_b2),
            _c(c),
            _c2(c*c),
            _cm2(1/_c2) {
            // ORSA_DEBUG("allocating (*this) = %x",this);
            _init();
        }
    private:
        void _init();
    protected:
        ~EllipsoidShape() {
            // ORSA_DEBUG("freeing (*this) = %x",this);
            gsl_root_fdfsolver_free (s);
        }
    
    public:
        ShapeType getType() const {
            return SHAPE_ELLIPSOID;
        }
    
    public:
        void getABC(double & a,
                    double & b,
                    double & c) const {
            a = _a;
            b = _b;
            c = _c;
        }

    public:
        double volume() const {
            return 4.0*orsa::pi()*_a*_b*_c/3.0;
        }
        
    public:
        bool isInside(const Vector &) const;
    
    public:
        const Vector closestVertex(const Vector &) const;
    
    protected:
        bool _updateCache() const;
    
    public:
        bool rayIntersection(orsa::Vector & intersectionPoint,
                             orsa::Vector & normal,
                             const orsa::Vector & P,
                             const orsa::Vector & u,
                             const bool fullLine = false) const; 
    
    public:
        // v should belong to the Ellipsoid, no test is performed
        inline orsa::Vector normalVector(const orsa::Vector v) const {
            return (orsa::Vector(v.getX()*_am2,
                                 v.getY()*_bm2,
                                 v.getZ()*_cm2).normalize());
        }
        //
        /* 
           inline orsa::Vector normalVector(const orsa::Vector v) const {
           orsa::Vector n(v.getX()*_am2,
           v.getY()*_bm2,
           v.getZ()*_cm2);
           n.normalize();
           return (n);
           }
        */
        
    protected:
        const double _a, _a2, _am2, _b, _b2, _bm2, _c, _c2, _cm2;
    public:
        mutable double closestVertexEpsilonRelative;
        // mutable double cV_l;
    public:
        struct cV_par {
            double a2, b2, c2;
            double Px, Py, Pz;
        };
    protected:
        mutable cV_par params;
        mutable gsl_function_fdf cV_FDF;
        const gsl_root_fdfsolver_type * T;
        gsl_root_fdfsolver * s;
        // double cV_S (const double & P, const double & l, const double & r2) const;
        // double cV_dS2_dl_over_r2 (const double & P, const double & l, const double & r2) const;
        // double cV_f (double l, void * params) const;
        // double cV_df (double l, void * params) const;
        // void cV_fdf (double l, void *params, double *y, double *dy) const;
    };
    
    //! Utility function: returns true if a half-line intersects a triangle
    //! Line starts at point P and with direction u 
    //! Triangle t1, t2, t3
    //! Use the bool to set whether the sign of u counts or not.
    bool rayIntersectsTriangle(orsa::Vector & intersectionPoint,
                               const orsa::Vector & P,
                               const orsa::Vector & u,
                               const orsa::Vector & t1,
                               const orsa::Vector & t2,
                               const orsa::Vector & t3,
                               const bool fullLine = false);
  
}; // namespace orsa

#endif // _ORSA_SHAPE_
