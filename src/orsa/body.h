#ifndef _ORSA_BODY_
#define _ORSA_BODY_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/cache.h>
#include <orsa/datetime.h>
#include <orsa/debug.h>
#include <orsa/massCluster.h>
#include <orsa/paulMoment.h>
#include <orsa/propulsion.h>
#include <orsa/quaternion.h>
// #include <orsa/reference_frame.h>
#include <orsa/vector.h>
#include <orsa/print.h>
#include <orsa/shape.h>

#include <string>

namespace orsa {
  
    class Body;
    class BodyInitialConditions;
  
    /***/
  
    class BodyProperty : public osg::Referenced {
    public:
        BodyProperty() : osg::Referenced(true) { }
    protected:
        virtual ~BodyProperty() { }
    public:
        //! BP_CONSTANT is not a function of time.
        //! BP_PRECOMPUTED is a funcion of time, but is already computed, useful for callbacks.
        //! BP_DYNAMIC is a function of time, and needs to be computed.
        enum BodyPropertyType {
            BP_CONSTANT,   
            BP_PRECOMPUTED, 
            BP_DYNAMIC   
        };
    public:
        virtual BodyPropertyType type() const = 0;
    public:
        virtual bool dynamic() const { return (type() == BP_DYNAMIC); }
    public:	
        /* update() should be implemented by all properties of type BP_PRECOMPUTED */
        virtual bool update(const orsa::Time &) = 0;
    public:
        // encapsulate the update and data retrevial in a lock/unlock, useful for properties of type BP_PRECOMPUTED in multi-threaded applications
        virtual void lock() = 0;
        virtual void unlock() = 0;
    };
    
    /***/
  
    /*
     * Three reference systems:
     * - global: the global, inertial frame; all coordinates are referred to global unless states otherwise
     * -  local: the frame of a body, where inertia matrix and PaulMoment are defined, and where the center of mass is always in (0,0,0)
     * -  shape: the frame where the shape model is defined, usually offset and rotated with respect to "local"
     *
     * To transform between "local" and "shape", a translation and a rotation are needed in general
     * rotation given by shapeToLocal and localToShape matrices
     * translation by the centerOfMass vector
     *
     *
     */
  
    /* 
     * centerOfMass is in the coordinates used for the shape
     * shapeToLocal is the rotation matrix from the shape system to the diagonal system ( = "local")
     * inertialMatrix is diagonal if shapeToLocal and localToShape are correct
     * paulMoment is about the centerOfMass and in the principalAxis system
     *
     * Usually, to compute these vars, a mass distribution is assumed,
     * or the data is otherwise obtained experimentally
     * 
     */
  
    class LocalShapeData : public osg::Referenced {
    public:
        LocalShapeData() : osg::Referenced(true) { }
    public:
        LocalShapeData(const LocalShapeData & data) : osg::Referenced(true) {
            localShape    = data.localShape;
            originalShape = data.originalShape;
            cm            = data.cm;
            s2l           = data.s2l;
        }
    public:
        LocalShapeData & operator = (const LocalShapeData & data) {
            localShape    = data.localShape;
            originalShape = data.originalShape;
            cm            = data.cm;
            s2l           = data.s2l;
            return (*this);
        }
    protected:
        virtual ~LocalShapeData() { } 
    public: 
        osg::ref_ptr<const orsa::Shape> localShape;
    public:
        osg::ref_ptr<const orsa::Shape> originalShape;
        orsa::Vector cm;
        orsa::Matrix s2l;
    };
  
    class LocalShapeCache : public osg::Referenced {
    public:
        LocalShapeCache() : osg::Referenced(true) { }
    public:
        LocalShapeCache(const LocalShapeCache & cache) : osg::Referenced(true) {
            data = cache.data;
        }
    public:
        LocalShapeCache & operator = (const LocalShapeCache & cache) {
            data = cache.data;
            return (*this);
        }
    protected:
        virtual ~LocalShapeCache() { }
    public:
        osg::ref_ptr<LocalShapeData> data;
    };
  
    class InertialBodyProperty : public BodyProperty {
    public:
        InertialBodyProperty() : BodyProperty() {
            // ORSA_DEBUG("creating object: %x",this);
            // IMPORTANT, each IPB needs one od this, to pass cached computations along
            ls_cache = new LocalShapeCache;
        }
    public:
        InertialBodyProperty(const InertialBodyProperty & ibp) : BodyProperty(ibp) {
            // ORSA_DEBUG("creating object: %x copy of %x",this,&ibp);
            ls_cache = ibp.ls_cache;
            // ORSA_DEBUG("ls_cache: %x",ls_cache.get());
        }
    public:
        InertialBodyProperty & operator = (const InertialBodyProperty & ibp) {
            ls_cache = ibp.ls_cache;
            // ORSA_DEBUG("ls_cache: %x",ls_cache.get());
            return (*this);
        }
    protected:
        virtual ~InertialBodyProperty() {
            // ORSA_DEBUG("destroying object: %x",this);
        }
    public:
        /* public: 
           InertialBodyProperty() : BodyProperty() {
           // ORSA_DEBUG("EMPTY CONSTRUCTOR");
           }
           public:
           InertialBodyProperty(const InertialBodyProperty & ibp) : BodyProperty() {
           // ORSA_DEBUG("COPY CONSTRUCTOR");
           // need to write explicitly this constructor to copy local members
           local_shape = ibp.local_shape;
           ls_prev_s   = ibp.ls_prev_s;
           ls_prev_cm  = ibp.ls_prev_cm;
           ls_prev_s2l = ibp.ls_prev_s2l;
           }
        */
    public:
        virtual double mass() const = 0;
        virtual const orsa::Shape * originalShape() const = 0;
        virtual orsa::Vector centerOfMass() const = 0;
        virtual orsa::Matrix shapeToLocal() const = 0;
        virtual orsa::Matrix localToShape() const = 0;
        virtual orsa::Matrix inertiaMatrix() const = 0;
        virtual const orsa::PaulMoment * paulMoment() const = 0;
        virtual const orsa::MassCluster * massCluster() const = 0;
    public:
        virtual bool setMass(const double &) = 0;
        virtual bool setOriginalShape(const orsa::Shape *) = 0;
        virtual bool setCenterOfMass(const orsa::Vector &) = 0;
        virtual bool setShapeToLocal(const orsa::Matrix &) = 0;
        virtual bool setLocalToShape(const orsa::Matrix &) = 0;
        virtual bool setInertiaMatrix(const orsa::Matrix &) = 0;
        virtual bool setPaulMoment(const orsa::PaulMoment *) = 0;
        virtual bool setMassCluster(const orsa::MassCluster *) = 0;
    public:
        // the shape, transformed into "local" coordinates
        virtual const orsa::Shape * localShape() const;
    protected:
        mutable osg::ref_ptr<LocalShapeCache> ls_cache;
    public:
        virtual InertialBodyProperty * clone() const = 0;
    };
  
    // NOTE: in general "operator =" of non-const IBP must copy "ls_cache" along. 
  
    class PointLikeConstantInertialBodyProperty : public InertialBodyProperty {
    public:
        PointLikeConstantInertialBodyProperty(const double & m) : 
            InertialBodyProperty(),
            _m(m) { }
    public:
        PointLikeConstantInertialBodyProperty(const PointLikeConstantInertialBodyProperty & ibp) : 
            InertialBodyProperty(ibp),
            _m(ibp._m) { }
    public:
        PointLikeConstantInertialBodyProperty & operator = (const PointLikeConstantInertialBodyProperty &) {
            ORSA_ERROR("this class is not supposed to change...");
            return (*this);
        }
    protected:
        const double _m;
    protected:
        double mass() const { return _m; }
        const orsa::Shape * originalShape() const { return 0; }
        orsa::Vector centerOfMass() const { return orsa::Vector(0,0,0); }
        orsa::Matrix shapeToLocal() const { return orsa::Matrix::identity(); }
        orsa::Matrix localToShape() const { return orsa::Matrix::identity(); }
        orsa::Matrix inertiaMatrix() const { return orsa::Matrix::identity(); }
        const orsa::PaulMoment * paulMoment() const { return 0; }
        const orsa::MassCluster * massCluster() const { return 0; }
    public:
        bool setMass(const double &) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
        bool setOriginalShape(const orsa::Shape *) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
        bool setCenterOfMass(const orsa::Vector &) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
        bool setShapeToLocal(const orsa::Matrix &) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
        bool setLocalToShape(const orsa::Matrix &) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
        bool setInertiaMatrix(const orsa::Matrix &) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
        bool setPaulMoment(const orsa::PaulMoment *) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
        bool setMassCluster(const orsa::MassCluster *) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
    public:
        InertialBodyProperty * clone() const {
            return new PointLikeConstantInertialBodyProperty(*this);
        }
    public:
        BodyPropertyType type() const { return BP_CONSTANT; }
    public:
        bool update(const orsa::Time &) { return true; }
    public:
        void lock() { }
        void unlock() { }
    };
    
    class ConstantInertialBodyProperty : public InertialBodyProperty {
    public:
        ConstantInertialBodyProperty(const double           & m,
                                     const orsa::Shape      * s,
                                     const orsa::Vector     & cm,
                                     const orsa::Matrix     & s2l,
                                     const orsa::Matrix     & l2s,
                                     const orsa::Matrix     & I,
                                     const orsa::PaulMoment * pm,
                                     const orsa::MassCluster * mc) :
            InertialBodyProperty(), 
            _m(m),
            _s(s),
            _cm(cm),
            _s2l(s2l),
            _l2s(l2s),
            _I(I),
            _pm(pm),
            _mc(mc) { }
    public:
        ConstantInertialBodyProperty(const ConstantInertialBodyProperty & ibp) :
            InertialBodyProperty(ibp), 
            _m(ibp._m),
            _s(ibp._s),
            _cm(ibp._cm),
            _s2l(ibp._s2l),
            _l2s(ibp._l2s),
            _I(ibp._I),
            _pm(ibp._pm),
            _mc(ibp._mc) { }
    public:
        ConstantInertialBodyProperty & operator = (const ConstantInertialBodyProperty &) {
            ORSA_ERROR("this class is not supposed to change...");
            return (*this);
        }
    protected:
        const double _m;
        osg::ref_ptr<const orsa::Shape> _s;
        const orsa::Vector _cm;
        const orsa::Matrix _s2l;
        const orsa::Matrix _l2s;
        const orsa::Matrix _I;
        osg::ref_ptr<const orsa::PaulMoment> _pm;
        osg::ref_ptr<const orsa::MassCluster> _mc;
    public:
        double mass() const { return _m; }
        const orsa::Shape * originalShape() const { return _s.get(); }
        orsa::Vector centerOfMass() const { return _cm; }
        orsa::Matrix shapeToLocal() const { return _s2l; }
        orsa::Matrix localToShape() const { return _l2s; }
        orsa::Matrix inertiaMatrix() const { return _I; }
        const orsa::PaulMoment * paulMoment() const { return _pm.get(); }
        const orsa::MassCluster * massCluster() const { return _mc.get(); }
    public:
        bool setMass(const double &) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
        bool setOriginalShape(const orsa::Shape *) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
        bool setCenterOfMass(const orsa::Vector &) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
        bool setShapeToLocal(const orsa::Matrix &) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
        bool setLocalToShape(const orsa::Matrix &) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
        bool setInertiaMatrix(const orsa::Matrix &) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
        bool setPaulMoment(const orsa::PaulMoment *) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
        bool setMassCluster(const orsa::MassCluster *) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
    public:
        InertialBodyProperty * clone() const {
            return new ConstantInertialBodyProperty(*this);
        }
    public:
        BodyPropertyType type() const { return BP_CONSTANT; }
    public:
        bool update(const orsa::Time &) { return true; }
    public:
        void lock() { }
        void unlock() { }
    };
    
    /***/
  
    class TranslationalBodyProperty : public BodyProperty {
    public:
        virtual orsa::Vector position() const = 0;
        virtual orsa::Vector velocity() const = 0;
    public:
        virtual bool setPosition(const orsa::Vector &) = 0;
        virtual bool setVelocity(const orsa::Vector &) = 0;
    public:
        virtual TranslationalBodyProperty * clone() const = 0;
    };
  
    class ConstantTranslationalBodyProperty : public TranslationalBodyProperty {
    public:
        ConstantTranslationalBodyProperty(const orsa::Vector & r) :
            TranslationalBodyProperty(),
            _position(r) { }
    public:
        BodyPropertyType type() const { return BP_CONSTANT; }
    public:
        bool update(const orsa::Time &) { return true; }
    public:
        void lock() { }
        void unlock() { }
    public:
        orsa::Vector position() const { 
            return _position;
        }
    public:
        orsa::Vector velocity() const {
            return orsa::Vector(0,0,0);
        }
   public:
        bool setPosition(const orsa::Vector &) { 
            ORSA_DEBUG("cannot change position");
            return false;
        }
    public: 
        bool setVelocity(const orsa::Vector &) {
            ORSA_DEBUG("cannot change velocity");
            return false;
        }
    public:
        TranslationalBodyProperty * clone() const {
            return new ConstantTranslationalBodyProperty(*this);
        }
    protected:
        const orsa::Vector _position;
    };
    
    class PrecomputedTranslationalBodyProperty : public TranslationalBodyProperty {
    public:
        BodyPropertyType type() const { return BP_PRECOMPUTED; }
    public:
        bool setPosition(const orsa::Vector &) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
    public:
        bool setVelocity(const orsa::Vector &) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
    public:
        TranslationalBodyProperty * clone() const {
            return (TranslationalBodyProperty *) this;
        }
    };
  
    class DynamicTranslationalBodyProperty : public TranslationalBodyProperty {
    public:
        BodyPropertyType type() const { return BP_DYNAMIC; }
    public:
        bool update(const orsa::Time &) { return true; }
    public:
        void lock() { }
        void unlock() { }
    public:
        orsa::Vector position() const { 
            // ORSA_DEBUG("_position: %.20e",_position.length());
            // orsa::print(_position);
            return _position;
        }
    public:
        orsa::Vector velocity() const {
            // ORSA_DEBUG("_velocity: %.20e",_velocity.length());
            // orsa::print(_velocity);
            return _velocity;
        }
    public:
        bool setPosition(const orsa::Vector & r) { 
            _position = r;
            // ORSA_DEBUG("_position: %.20e",_position.length());
            // orsa::print(_position);
            return true;
        }
    public: 
        bool setVelocity(const orsa::Vector & v) {
            _velocity = v;
            // ORSA_DEBUG("_velocity: %.20e",_velocity.length());
            // orsa::print(_velocity);
            return true;
        }
    public:
        TranslationalBodyProperty * clone() const {
            return new DynamicTranslationalBodyProperty(*this);
        }
    protected:
        orsa::Cache<orsa::Vector> _position;
        orsa::Cache<orsa::Vector> _velocity;
    };
  
    /***/
  
    class RotationalBodyProperty : public BodyProperty {
    public:
        virtual bool get(orsa::Quaternion & q,
                         orsa::Vector     & omega) const = 0;
    public:
        virtual orsa::Quaternion getQ()    const = 0;
        virtual orsa::Vector     getOmega() const = 0;
    
    public:
        /* 
           orsa::Quaternion qRotated(const orsa::Vector & omega,
           const orsa::Time   & dt) const;
        */
    
    public:
        virtual bool set(const orsa::Quaternion & q,
                         const orsa::Vector     & omega) = 0;
    public:
        virtual RotationalBodyProperty * clone() const = 0;
  
    public:
        static orsa::Quaternion qDot (const orsa::Quaternion & q,
                                      const orsa::Vector     & omega);
    
        static orsa::Quaternion qDotDot (const orsa::Quaternion & q,
                                         const orsa::Vector     & omega,
                                         const orsa::Vector     & omegaDot);
    
        static orsa::Vector omega (const orsa::Quaternion & q,
                                   const orsa::Quaternion & qDot);
    
        static orsa::Vector omegaDot (const orsa::Quaternion & q,
                                      const orsa::Quaternion & qDot,
                                      const orsa::Quaternion & qDotDot);
    
        static orsa::Quaternion qFiniteRotation (const orsa::Quaternion & q,
                                                 const orsa::Vector     & omega,
                                                 const orsa::Time       & dt);
    
        static orsa::Vector newOmega (const orsa::Vector & omega,
                                      const orsa::Vector & omegaDot,
                                      const orsa::Time   & dt);
    };
  
    class ConstantRotationalBodyProperty : public RotationalBodyProperty {
    public:
        BodyPropertyType type() const { return BP_CONSTANT; }
    public:
        bool update(const orsa::Time &) { return true; }
    public:
        void lock() { }
        void unlock() { }
    };
    
    class PrecomputedRotationalBodyProperty : public RotationalBodyProperty {
    public:
        BodyPropertyType type() const { return BP_PRECOMPUTED; }
    public:
        virtual bool set(const orsa::Quaternion &,
                         const orsa::Vector     &) {
            ORSA_ERROR("this method should not have been called, please check your code.");
            return false;
        }
    public:
        RotationalBodyProperty * clone() const {
            return (RotationalBodyProperty *) this;
        }
    };
  
    class DynamicRotationalBodyProperty : public RotationalBodyProperty {
    public:
        BodyPropertyType type() const { return BP_DYNAMIC; }
    public:
        bool update(const orsa::Time &) { return true; }
    public:
        void lock() { }
        void unlock() { }
    public:
        bool get(orsa::Quaternion & q,
                 orsa::Vector     & omega) const {
            q     = _q;
            omega = _omega;
            return true;
        }
    public:
        orsa::Quaternion getQ()     const { return _q; }
        orsa::Vector     getOmega() const { return _omega; }
    public:
        bool set(const orsa::Quaternion & q,
                 const orsa::Vector     & omega) {
            _q     = q;
            _omega = omega;
            // ORSA_DEBUG("omega: %Fg",_omega.length());
            return true;
        }
    public:
        RotationalBodyProperty * clone() const {
            return new DynamicRotationalBodyProperty(*this);
        }
    protected:
        orsa::Cache<orsa::Quaternion> _q;
        orsa::Cache<orsa::Vector>     _omega;
    };
  
    /***/
    
    class UpdateIBPS : public osg::Referenced {
    public:
        UpdateIBPS() : osg::Referenced(true) { }
    protected:
        virtual ~UpdateIBPS() { }
    public:
        virtual bool update(const orsa::Time & t,
                            InertialBodyProperty * inertial,
                            TranslationalBodyProperty * translational,
                            RotationalBodyProperty * rotational) {
            if (inertial) {
                inertial->update(t);
            }
            if (translational) {
                translational->update(t);
            }
            if (rotational) {
                rotational->update(t);
            }
            return true;
        }
    public:
        virtual UpdateIBPS * clone() const {
            return new UpdateIBPS(*this);
        }
    };
    
    /***/
    
    class IBPS {
    public:
        IBPS();
    public:
        IBPS(const IBPS &);
    public:
        virtual ~IBPS();
    public:
        orsa::Cache<orsa::Time> time;
    public:
        osg::ref_ptr<InertialBodyProperty> inertial;
    public:
        osg::ref_ptr<TranslationalBodyProperty> translational;
    public:
        osg::ref_ptr<RotationalBodyProperty> rotational;
    public:
        osg::ref_ptr<UpdateIBPS> updateIBPS;
    public:
        IBPS & operator = (const IBPS &);
        
    public:
        virtual bool update(const orsa::Time & t) {
            return updateIBPS->update(t,inertial,translational,rotational);
        }
        
    public:
        void lock() {
            if (inertial.get()) {
                inertial->lock();
            }
            if (translational.get()) {
                translational->lock();
            }
            if (rotational.get()) {
                rotational->lock();
            }
        }
    public:
        void unlock() {
            if (inertial.get()) {
                inertial->unlock();
            }
            if (translational.get()) {
                translational->unlock();
            }
            if (rotational.get()) {
                rotational->unlock();
            }
        }
        
    public:
        bool dynamic() const {
            if (inertial.get()) {
                if (inertial->dynamic()) return true;
            }
            if (translational.get()) {
                if (translational->dynamic()) return true;
            }
            if (rotational.get()) {
                if (rotational->dynamic()) return true;
            }
            return false;
        }
    
    public:
        bool tmp;
    
        // time ordering, needed when used in orsa::Interval<>
    public:
        inline bool operator == (const IBPS & rhs) const {
            return (time == rhs.time);
        }
    public:
        inline bool operator != (const IBPS & rhs) const {
            return (time != rhs.time);
        }
    public:
        inline bool operator < (const IBPS & rhs) const {
            return (time < rhs.time);
        }
    public:
        inline bool operator > (const IBPS & rhs) const {
            return (time > rhs.time);
        }
    public:
        inline bool operator <= (const IBPS & rhs) const {
            return (time <= rhs.time);
        }
    public:
        inline bool operator >= (const IBPS & rhs) const {
            return (time >= rhs.time);
        }
    };
  
    typedef orsa::IBPS InstantaneousBodyPropertySet;
  
    /***/
  
    class Body : public osg::Referenced {
    public:	
        Body();
    private:
        void _init();
    protected:
        virtual ~Body();
    
    public:
        typedef unsigned int BodyID;
        BodyID id() const { return _id; }
    private:
        static BodyID bodyID_counter;
        const BodyID _id;
        
    public:
        bool setInitialConditions(const orsa::IBPS & ibps) {
            // these checks do not completely enforce presence of inertial and translational components, but help indeed achieve this goal.
            if (ibps.inertial.get() == 0) {
                ORSA_DEBUG("ibps.inertial must be set.");
                return false;
            }	 
            if (ibps.translational.get() == 0) {
                ORSA_DEBUG("ibps.translational must be set.");
                return false;
            }	
            _ibps = ibps;
            return true;
        }
        
    public:
        const orsa::IBPS & getInitialConditions() const {
            return _ibps;
        }
    protected:
        orsa::Cache<orsa::IBPS> _ibps;
        
    public:
        virtual bool setName(const std::string &);
        virtual const std::string & getName() const;
    protected:
        std::string _name;
        
    public:
        osg::ref_ptr<const orsa::Propulsion> propulsion;
    
        // comet particle beta
    public:
        mutable orsa::Cache<double> beta;
        osg::ref_ptr<const orsa::Body> betaSun;
        
        // bio
    public:
        mutable orsa::Cache<orsa::Time> birthTime;
        mutable orsa::Cache<orsa::Time> deathTime;
    public:
        bool alive(const orsa::Time & t) const {
            if (birthTime.isSet()) {
                if (t < birthTime) return false;
            }
            if (deathTime.isSet()) {
                if (t > deathTime) return false;
            }
            return true;
        }
    
        // light source
        // should add a 'graphics' entry to IBPS, for this and other related entries
    public:
        orsa::Cache<bool> isLightSource;
    
    public:
        //! all the bodies with the same value of nonInteractingGroupID do not interact with each other
        //! example: particles with small mass, that don't see each other, but that all interact with sun and planets
        // orsa::Cache<unsigned int> nonInteractingGroupID;
        //! simpler (and faster): one single group for all non-interacting particles (cannot think of any case where more than one group is really needed)
        //! defaults to false: iteracts with all other bodies
        bool nonInteractingGroup;
    
    };
  
} // namespace orsa

#endif // _ORSA_BODY_
