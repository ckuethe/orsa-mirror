#ifndef _ORSA_INTERACTION_
#define _ORSA_INTERACTION_

#include <orsa/datetime.h>
#include <orsa/vector.h>
#include <osg/Referenced>
#include <osg/ref_ptr>
#include <omp.h>

namespace orsa {

    class Body;
    class BodyGroup;
    class BodyPair;
    class PaulMoment;

    class Interaction : public osg::Referenced {

    public:
        Interaction();
    protected:
        virtual ~Interaction() { }

    public:
        typedef std::vector<orsa::Vector> InteractionVector;

    public:
        bool accelerationAndTorque(InteractionVector & a,
                                   InteractionVector & N,
                                   orsa::BodyGroup   * bg,
                                   const orsa::Time  & t,
                                   bool computeAcceleration = true,
                                   bool computeTorque = true) const;
        /* bool acceleration(InteractionVector & a,
	   orsa::BodyGroup   * bg,
	   const orsa::Time  & t) const;
	   bool torque(InteractionVector & N,
	   orsa::BodyGroup   * bg,
	   const orsa::Time  & t) const;
	*/
	
    public:
        bool dependsOnVelocity() const { return true; }

    protected:
        bool bodyPairAccelerationTerm(orsa::Vector     & a_ref_b,
                                      orsa::Vector     & a_b,
                                      orsa::BodyGroup  * bg,
                                      const BodyPair   * bp,
                                      const orsa::Time & t) const;

    protected:
        osg::ref_ptr<PaulMoment> dummyPaulMoment;

    };

} // namespace orsa

#endif // _ORSA_INTERACTION_





//#ifndef _ORSA_INTERACTION_
//#define _ORSA_INTERACTION_
//
//#include <orsa/datetime.h>
//#include <orsa/vector.h>
//
//#include <osg/Referenced>
//#include <osg/ref_ptr>
//
//namespace orsa {
//
//    class Body;
//    class BodyGroup;
//    class BodyPair;
//    class PaulMoment;
//
//    class Interaction : public osg::Referenced {
//
//        /*
//           public:
//           enum InteractionType {
//           NEWTON = 1,
//           FAST_RELATIVISTIC_CORRECTIONS = 2,
//           RELATIVISTIC_CORRECTIONS = 4,
//           MULTIPOLES = 8,
//           TREE = 16,
//           DEFAULT = NEWTON
//           };
//        */
//
//    public:
//        Interaction();
//    protected:
//        virtual ~Interaction() { }
//
//    public:
//        // #warning "SHOULD rename it, no longer BodyID aligned, but bg::BodyList index"
//        // typedef std::vector<orsa::Vector> BodyIDVector;
//        typedef std::vector<orsa::Vector> InteractionVector;
//    public:
//        bool acceleration(InteractionVector & a,
//                          orsa::BodyGroup   * bg,
//                          const orsa::Time  & t) const;
//    public:
//        // this should be protected:
//        bool bodyPairAccelerationTerm(orsa::Vector     & a_ref_b,
//                                      orsa::Vector     & a_b,
//                                      orsa::BodyGroup  * bg,
//                                      const BodyPair   * bp,
//                                      const orsa::Time & t) const;
//
//    public:
//        bool torque(InteractionVector & N,
//                    orsa::BodyGroup   * bg,
//                    const orsa::Time  & t) const;
//
//        /*
//           public:
//           virtual bool setType(const InteractionType) = 0;
//           virtual InteractionType getType() const = 0;
//        */
//
//    public:
//        // virtual bool dependsOnVelocity() const = 0;
//        bool dependsOnVelocity() const { return true; }
//
//    protected:
//        osg::ref_ptr<PaulMoment> dummyPaulMoment;
//
//    };
//
//} // namespace orsa
//
//#endif // _ORSA_INTERACTION_
