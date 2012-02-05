#include <orsa/interaction.h>
#include <orsa/bodygroup.h>
#include <orsa/paul.h>
#include <orsa/print.h>
#include <orsa/util.h>

using namespace orsa;

// Constructor
Interaction::Interaction() : osg::Referenced(true) {
    dummyPaulMoment = new PaulMoment(0);
    dummyPaulMoment->setM(1,0,0,0);
}

// Acceleration (A) and Torque (T)
bool Interaction::accelerationAndTorque(InteractionVector & a,
                                        InteractionVector & N,
                                        orsa::BodyGroup   * bg,
                                        const orsa::Time  & t,
                                        bool computeAcceleration,
                                        bool computeTorque) const {
    
    // Something to do for A or T
    if (computeAcceleration || computeTorque) {
        
        // Initialization A & T
        BodyGroup::BodyList bl = bg->getBodyList();
        unsigned int blSize = bl.size();
        {
            if (computeAcceleration && computeTorque) {
                if (a.size() != blSize) a.resize(blSize);
                if (N.size() != blSize) N.resize(blSize);
                for (unsigned int k=0; k<blSize; ++k) {
                    a[k] = orsa::Vector(0,0,0);
                    N[k] = orsa::Vector(0,0,0);
                }
            } else if (!computeTorque) {
                if (a.size() != blSize) a.resize(blSize);
                for (unsigned int k=0; k<blSize; ++k) {
                    a[k] = orsa::Vector(0,0,0);
                }                
            } else if (!computeAcceleration) {
                if (N.size() != blSize) N.resize(blSize);
                for (unsigned int k=0; k<blSize; ++k) {
                    N[k] = orsa::Vector(0,0,0);
                }
            }
        }
        
        // Loop part for A & T
        {
            // Loop global variables
            double m_ref_b, m_b;            // A & T
            
            // Loop local variables
            IBPS ref_b_ibps;                // A & T
            IBPS b_ibps;                    // A & T
            orsa::Matrix ref_b_g2l;         // A & T
            orsa::Matrix b_g2l;             // A & T
            orsa::Vector R;                 // A & T
            orsa::Vector thrust;            // A
            double b_radius;                // A
            double ref_b_radius;            // A
            orsa::Vector _d;                // A
            double _l;                      // A
            orsa::Vector accTerm;           // A 
            bool compute_ajk;               // A
            orsa::Vector torqueTerm;        // T
            
            bool continueA;                 // A
            bool continueT;                 // T
            
            // Loop on bodies for A & T
#pragma omp parallel for private(continueA,continueT,m_ref_b,m_b,ref_b_ibps,b_ibps,ref_b_g2l,R,thrust,b_radius,ref_b_radius,_d,_l,accTerm,compute_ajk,torqueTerm) schedule(dynamic)
            for (unsigned int j=0; j<blSize; ++j) {
	      
                // Declaration
                const orsa::Body * ref_b = bl[j].get();
                continueA = false;
                continueT = false;
                
                // Current body initialization
                {

                    // A & T
                    {
                        if (!ref_b->alive(t)) {
                            std::cout<<"alive"<<std::endl;
                            continue;
                        }
                    } // A & T

                    // A 
                    {
                        if (ref_b->getInitialConditions().translational.get()) {
                            if (!(ref_b->getInitialConditions().translational->dynamic())) {
                                continueA = true;
                            }
                        } else {
                            continueA = true;
                        }
                    } // A
                    
                    // T
                    {
                        if (ref_b->getInitialConditions().rotational.get()) {
                            if (!(ref_b->getInitialConditions().rotational->dynamic())) {
                                continueT = true;
                            }
                        } else {
                            continueT = true;
                        }
                    }
                    
                    // A & T   
                    {
                        if (continueA && continueT) {
                            continue;
                        }
                        
                        if (!(bg->getInterpolatedIBPS(ref_b_ibps,
                                                      ref_b,
                                                      t))) {
                            ORSA_DEBUG("problems...");
                            //return false;
                        }
                        
                        if (!bg->getInterpolatedMass(m_ref_b,ref_b,t)) {
                            ORSA_DEBUG("problems...");
                        }
                    } // A & T
                    
                    // A
                    if (computeAcceleration && !continueA) {
                        thrust.set(0,0,0);
                        if (ref_b->propulsion.get()) {
                            if (!dependsOnVelocity()) {
                                ORSA_DEBUG("PROBLEMS: when using a Propulsion, in general you need to set orsa::Interaction::dependsOnVelocity() as TRUE");
                            }
                            thrust = ref_b->propulsion->getThrust(t);
                        }
                    } // A
                    
                    // T
                    if (computeTorque && !continueT) {
                        osg::ref_ptr<const PaulMoment> ref_b_pm =
                        (ref_b_ibps.inertial->paulMoment()) ?
                        (ref_b_ibps.inertial->paulMoment()) :
                        (dummyPaulMoment.get());
                    } // T
               
                } // Current body initialization  
                
                // Internal loop on other bodies
                for (unsigned int k=0; k<j; ++k) {
                   
                    // Declaration
                    const orsa::Body * b = bl[k].get();
                    
                    // Internal loop initialization
                    {
                        // A & T
                        {
                            compute_ajk = false;

                            if (!b->alive(t)) {
                                continue;
                            }

                            if (!bg->getInterpolatedMass(m_b,b,t)) {
                                ORSA_DEBUG("problems...");
                            }

                            if (ref_b->nonInteractingGroup &&
                                b->nonInteractingGroup) {
                                continue;
                            }

                            if ((m_b == 0) &&
                                (m_ref_b == 0)) {
                                continue;
                            }

                            if (!(bg->getInterpolatedIBPS(b_ibps,
                                                          b,
                                                          t))) {
                                ORSA_DEBUG("problems...");
                                //return false;
                            }
                        } // A & T
                    } // Internal loop initialization
                    
                    // Main part of internal loop
                    {
                        if (ref_b_ibps.inertial->paulMoment() || b_ibps.inertial->paulMoment()) {

                            // Declaration
                            osg::ref_ptr<const PaulMoment> ref_b_pm =
                                    (ref_b_ibps.inertial->paulMoment()) ?
                                    (ref_b_ibps.inertial->paulMoment()) :
                                    (dummyPaulMoment.get());
                            osg::ref_ptr<const PaulMoment> b_pm =
                                    (b_ibps.inertial->paulMoment()) ?
                                    (b_ibps.inertial->paulMoment()) :
                                    (dummyPaulMoment.get());
                            
                            // A
                            if (computeAcceleration && !continueA) {

                                ref_b_g2l = orsa::globalToLocal(ref_b,bg,t);
                                b_g2l = orsa::globalToLocal(b,bg,t);

                                R =
                                    b_ibps.translational->position() -
                                    ref_b_ibps.translational->position();

                                b_radius =     b_ibps.inertial->localShape() ?     b_ibps.inertial->localShape()->boundingRadius() : 0;
                                ref_b_radius = ref_b_ibps.inertial->localShape() ? ref_b_ibps.inertial->localShape()->boundingRadius() : 0;

                                if ((b_radius+ref_b_radius) > R.length()) {

                                    ORSA_DEBUG("bodies too close: R<R1+R2, R=%g R1=%g R2=%g [km] [b1:%s] [b2:%s]",
                                               orsa::FromUnits(R.length(),orsa::Unit::KM,-1),
                                               orsa::FromUnits(b_radius,orsa::Unit::KM,-1),
                                               orsa::FromUnits(ref_b_radius,orsa::Unit::KM,-1),
                                               b->getName().c_str(),
                                               ref_b->getName().c_str());
                                    ORSA_DEBUG("reverting to pointlike...");
#warning should handle this better....

                                    _d =
                                        b_ibps.translational->position() -
                                        ref_b_ibps.translational->position();

                                    _l = _d.length();

                                    _d /= (_l*_l*_l);

                                    if (ref_b->betaSun == b) {
                                        accTerm = (1 - ref_b->beta) * _d;
                                        compute_ajk = true;
                                    } else {
                                        accTerm = _d;
                                        compute_ajk = true;
                                    }

                                } else {
                                    accTerm =
                                        Paul::gravitationalForce(ref_b_pm.get(),
                                                                 ref_b_g2l,
                                                                 b_pm.get(),
                                                                 b_g2l,
                                                                 R);

                                    compute_ajk = true;
                                }
                            } // A
                            
                            // T
                            if (computeTorque && !continueT) {
                                if ( (b_ibps.translational.get() &&
                                      ref_b_ibps.translational.get()) &&
                                     (b_ibps.rotational.get() &&
                                      ref_b_ibps.rotational.get()) ) {

                                    ref_b_g2l = orsa::globalToLocal(ref_b,bg,t);
                                    b_g2l = orsa::globalToLocal(b,bg,t);

                                    R =
                                        b_ibps.translational->position() -
                                        ref_b_ibps.translational->position();

                                    torqueTerm =
                                        Paul::gravitationalTorque(ref_b_pm.get(),
                                                                  ref_b_g2l,
                                                                  b_pm.get(),
                                                                  b_g2l,
                                                                  R);

                                    N[j] += orsa::Unit::G() * m_b     * torqueTerm;
                                    N[k] -= orsa::Unit::G() * m_ref_b * torqueTerm;

                                } else {
                                    if (!ref_b_ibps.translational.get()) ORSA_DEBUG("problems: [%s].translational.get() = %x",
                                                                                    ref_b->getName().c_str(),
                                                                                    ref_b_ibps.translational.get());
                                    if (!b_ibps.translational.get()) ORSA_DEBUG("problems: [%s].translational.get() = %x",
                                                                                b->getName().c_str(),
                                                                                b_ibps.translational.get());
                                }     
                            } // T
                            

                        } else {
                            
                            // A
                            if (computeAcceleration && !continueA) {
                                _d =
                                    b_ibps.translational->position() -
                                    ref_b_ibps.translational->position();

                                _l = _d.length();

                                if (_l > epsilon()) {

                                    _d /= (_l*_l*_l);

                                    if (ref_b->betaSun == b) {
                                        accTerm = (1 - ref_b->beta) * _d;
                                        compute_ajk = true;
                                    } else {
                                        accTerm = _d;
                                        compute_ajk = true;
                                    }

                                } else {

                                    ORSA_WARNING("skipping: zero distance between bodies [%s] and [%s]",
                                                 ref_b->getName().c_str(),
                                                 b->getName().c_str());
                                    ORSA_DEBUG("[%s] position:",b->getName().c_str());
                                    orsa::print(b_ibps.translational->position());
                                    ORSA_DEBUG("[%s] position:",ref_b->getName().c_str());
                                    orsa::print(ref_b_ibps.translational->position());
                                }
                            } // A
                            
                        }
                    } // Main part of internal loop
                
                    // Internal loop finalization
                    {
                        // A
                        if (computeAcceleration && !continueA) {
                            if (compute_ajk) {
                                a[j] += orsa::Unit::G() * m_b     * accTerm;
                                a[k] -= orsa::Unit::G() * m_ref_b * accTerm;
                            }
                        } // A
                    } // Internal loop finalization
                    
                } // Internal loop on other bodies

                // After internal loop
                {
                    // A
                    if (computeAcceleration && !continueA) {
                        if (ref_b->propulsion.get()) {
                            if (m_ref_b < orsa::epsilon()) {
                                ORSA_DEBUG("Propulsion problems: non-positive mass...");
                            } else {
                                a[j] += thrust / m_ref_b;
                            }
                        }
                    } // A
                } // After internal loop
                
            } // Loop on bodies for A & T 
        } // Loop part for A & T
    } // Something to do for A or T
    
    return true;
}

