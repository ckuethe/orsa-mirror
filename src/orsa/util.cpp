#include <orsa/util.h>
#include <orsa/print.h>
#include <orsa/statistic.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

using namespace orsa;

std::string & orsa::removeLeadingAndTrailingSpaces(std::string & s) {
    const int first = s.find_first_not_of(" ");
    s.erase(0,first);
    const int last  = s.find_last_not_of(" ");
    s.erase(last+1,s.size());
    return s;
}

std::string & orsa::removeAllSpaces(std::string & s) {
    s.erase(remove(s.begin(),s.end(),' '),s.end());
    return s;
}

std::string & orsa::removeLeadingPlusSign(std::string & s) {
    const int first = s.find_first_not_of("+");
    s.erase(0,first);
    return s;
}


/* This code was used in orsa/interaction.cpp at some point
   to determine on the fly a good value for dt. 
   The only unsolved problem is that of using non-tmp IBPS entries,
   as temporary entries are just wrong outside the integrator.
   
   if (ref_b->propulsion.get()) {
   const orsa::Vector propulsionThrust = ref_b->propulsion->getThrust(t);
   if (propulsionThrust.lengthSquared() > orsa::epsilon()*orsa::epsilon()) {
   const orsa::BodyGroup::BodyInterval * bi =
   bg->getBodyInterval(ref_b);
   const orsa::Time dt = 
   std::min(orsa::Time(0,0,0,5,0),
   std::min(bi->max().time-t,
   t-bi->min().time));
   // REMEMBER: now thrust is in FrenetSerret components
   orsa::Vector T, N, B;
   ORSA_DEBUG("bi->min().time.tmp: %i",bi->min().tmp);
   ORSA_DEBUG("bi->max().time.tmp: %i",bi->max().tmp);
   if (t > bi->min().time) {
   FrenetSerret(ref_b, bg,
   t,
   -dt,
   T, N, B);
   } else if (t < bi->max().time) {
   FrenetSerret(ref_b, bg,
   t,
   dt,
   T, N, B);
   } else {
   ORSA_DEBUG("interval smaller than dt");
   ORSA_DEBUG("--BODY-INTERVAL-TIME-RANGE--");
   print(bi->min().time);
   print(bi->max().time);
   ORSA_DEBUG("call time:");
   print(t);
   //
   T = orsa::Vector(1,0,0);
   N = orsa::Vector(0,1,0);
   B = orsa::Vector(0,0,1);
   }
   thrust = 
   propulsionThrust.getX()*T +	
   propulsionThrust.getY()*N +
   propulsionThrust.getZ()*B;
   }
   }
*/

bool orsa::FrenetSerret(const orsa::Body * b,
                        orsa::BodyGroup  * bg,
                        const orsa::Time & t,
                        const orsa::Time & dt,
                        orsa::Vector & T,
                        orsa::Vector & N,
                        orsa::Vector & B) {
    if (dt.getMuSec() == 0) return false;
    orsa::Vector r_t,   v_t;
    orsa::Vector r_tdt, v_tdt;
    if (bg->getInterpolatedPosVel(r_t,  v_t,  b,t) && 
        bg->getInterpolatedPosVel(r_tdt,v_tdt,b,t+dt)) {
        T = v_t.normalized();
        // const orsa::Vector dv = v_tdt-v_t;
        // should chech that velocity is not constant...
        N = ((v_tdt-v_t)/dt.get_d()).normalized();
        B = externalProduct(T,N).normalized();
        //
        /* 
           ORSA_DEBUG("b: [%s]",b->getName().c_str());
           print(t);
           print(dt);
           print(T);
           print(N);
           print(B);
        */
        //
        return true;
    } else {
        // ORSA_DEBUG("problems...");
        return false;
    }
}

bool orsa::eulerAnglesToMatrix(orsa::Matrix       & m,
                               const double & psi,
                               const double & theta,
                               const double & phi) {
  
    // ORSA_DEBUG("this code needs to be verified!!");
    // #warning "this code needs to be verified!!"
  
    m = orsa::Matrix::identity();
  
    m.rotZ(phi);
  
    m.rotX(theta);
  
    m.rotZ(psi);
  
    return true;
}	

bool orsa::matrixToEulerAngles(double       & psi,
                               double       & theta,
                               double       & phi,
                               const orsa::Matrix & m) {
  
    // const orsa::Matrix l2g = localToGlobal(t);
  
    const double cosTheta = m.getM33();
  
    // if (fabs(sinTheta) > epsilon()) {
    if ((1-cosTheta*cosTheta) > (epsilon()*epsilon())) {
    
        psi   = atan2(-m.getM13(),
                      m.getM23());
    
        phi   = atan2(-m.getM31(), 
                      -m.getM32());
    
        // should check this better...
        if (fabs(m.getM23()) > epsilon()) {
            theta = atan2(-m.getM23()/cos(psi),
                          m.getM33());
        } else if (fabs(m.getM31()) > epsilon()) {
            theta = atan2(m.getM31()/sin(phi),
                          m.getM33());
        } else {
            // we should not be here...
            ORSA_DEBUG("problems...");
            theta = 0;
        }
    
    } else {
    
        // ORSA_DEBUG("using singular code...");
    
        psi = theta = 0;
    
        phi = atan2(m.getM21(), 
                    m.getM11());
    
    }
  
    return true;
}	


orsa::Matrix orsa::QuaternionToMatrix (const orsa::Quaternion & q) {
  
    const double       s = q.getScalar();
    const orsa::Vector v = q.getVector();
  
    const double q0 = s;
    const double q1 = v.getX();
    const double q2 = v.getY();
    const double q3 = v.getZ();
  
    orsa::check(q0);
    orsa::check(q1);
    orsa::check(q2);
    orsa::check(q3);
  
    // Eq. (7.7) book by J.B. Kuipers on Quaternions
    return orsa::Matrix(2*q0*q0-1+2*q1*q1,   2*q1*q2+2*q0*q3, 2*q1*q3-2*q0*q2,
                        2*q1*q2-2*q0*q3,   2*q0*q0-1+2*q2*q2, 2*q2*q3+2*q0*q1,
                        2*q1*q3+2*q0*q2,     2*q2*q3-2*q0*q1, 2*q0*q0-1+2*q3*q3);
}

orsa::Quaternion orsa::MatrixToQuaternion (const orsa::Matrix & m) {
  
    if (fabs(m.determinant()-1) > 6*orsa::epsilon()) {
        ORSA_DEBUG("problem: call with a rotation matrix only: |det(m)-1| = %g",fabs(m.determinant()-1));
        return orsa::Quaternion();
    }
  
    // code discussion at http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
    // first, the magnitude
    double q0 = 0.5*sqrt(std::max(0.0,1+m.getM11()+m.getM22()+m.getM33()));
    double q1 = 0.5*sqrt(std::max(0.0,1+m.getM11()-m.getM22()-m.getM33()));
    double q2 = 0.5*sqrt(std::max(0.0,1-m.getM11()+m.getM22()-m.getM33()));
    double q3 = 0.5*sqrt(std::max(0.0,1-m.getM11()-m.getM22()+m.getM33()));
    // then the sign
    q1 = copysign(q1,m.getM23()-m.getM32());
    q2 = copysign(q2,m.getM31()-m.getM13());
    q3 = copysign(q3,m.getM12()-m.getM21());
  
    if (0) {
        // debug: compare new code with old one
        const double old_q1 = (m.getM23()-m.getM32())/(4*q0);
        const double old_q2 = (m.getM31()-m.getM13())/(4*q0);
        const double old_q3 = (m.getM12()-m.getM21())/(4*q0);
        ORSA_DEBUG("dq1/q1 = %+e   dq2/q2 = %+e   dq3/q3 = %+e",
                   (q1-old_q1)/fabs(q1),
                   (q2-old_q2)/fabs(q2),
                   (q3-old_q3)/fabs(q3));
    }
  
    return orsa::Quaternion(q0, orsa::Vector(q1,q2,q3));
}

orsa::Matrix orsa::localToGlobal(const orsa::Body       * b,
                                 const orsa::BodyGroup  * bg,
                                 const orsa::Time       & t) {
  
    if ((bg==0) || (b==0)) {
        ORSA_ERROR("invalid pointer...");
        return orsa::Matrix::identity();
    }
  
    orsa::IBPS ibps;
    if (bg->getInterpolatedIBPS(ibps, b, t)) { 
        if (ibps.rotational.get()) {
            ibps.rotational->lock();
            ibps.rotational->update(t);
            const orsa::Matrix m = QuaternionToMatrix(ibps.rotational.get()->getQ());
            ibps.rotational->unlock();
            return m;
        } else {
            // ORSA_DEBUG("problems... body: [%s]",b->getName().c_str());
            return orsa::Matrix::identity();
        }
    } else {
        ORSA_DEBUG("problems...body: [%s]",b->getName().c_str());
        return orsa::Matrix::identity();
    }
}

orsa::Matrix orsa::globalToLocal(const orsa::Body       * b,
                                 const orsa::BodyGroup  * bg,
                                 const orsa::Time       & t) {
  
    if ((bg==0) || (b==0)) {
        ORSA_ERROR("invalid pointer...");
        return orsa::Matrix::identity();
    }
  
    orsa::IBPS ibps;
    if (bg->getInterpolatedIBPS(ibps, b, t)) { 
        if (ibps.rotational.get()) {
            ibps.rotational->lock();
            ibps.rotational->update(t);
            const orsa::Matrix m = QuaternionToMatrix(ibps.rotational.get()->getQ());
            orsa::Matrix m_tr;
            orsa::Matrix::transpose(m, m_tr);
            ibps.rotational->unlock();
            return m_tr;
        } else {
            // ORSA_DEBUG("problems... body: [%s]",b->getName().c_str());
            return orsa::Matrix::identity();
        }
    } else {
        ORSA_DEBUG("problems... body: [%s]",b->getName().c_str());
        return orsa::Matrix::identity();
    }
}

// magnitude function
// alpha = solar phase angle = angle Sun-Asteroid-Observer
// G = slope parameter (G ~= 0.15)
double orsa::P (const double & alpha, 
                const double & G) {
    // ORSA_DEBUG("P:   alpha = %f",alpha.get_mpf_t());
    const double phi_1 = exp(-3.33*pow(tan(0.5*alpha),0.63));
    const double phi_2 = exp(-1.87*pow(tan(0.5*alpha),1.22));
    /* 
       ORSA_DEBUG("P = %f   alpha: %f   p1: %f   p2: %f",
       -2.5*log10((1.0-G)*phi_1+G*phi_2),
       alpha.get_mpf_t(),
       phi_1,
       phi_2);
    */
    return (-2.5*log10((1.0-G)*phi_1+G*phi_2));
}

double orsa::apparentMagnitude(const double & H,
                               const double & G,
                               const double & phaseAngle,
                               const double & neo2obs,
                               const double & neo2sun) {
    
    const double V = H + P(phaseAngle,G) + 
        5*log10(FromUnits(neo2obs,orsa::Unit::AU,-1)*FromUnits(neo2sun,orsa::Unit::AU,-1));
    
    return V;
}

double orsa::absoluteMagnitude(const double & V,
                               const double & G,
                               const double & phaseAngle,
                               const double & neo2obs,
                               const double & neo2sun) {
    
    const double H = V - P(phaseAngle,G) - 
        5*log10(FromUnits(neo2obs,orsa::Unit::AU,-1)*FromUnits(neo2sun,orsa::Unit::AU,-1));
    
    return H;
}

double orsa::asteroidDiameter(const double & p, 
                              const double & H) {
    if (p<0) {
        ORSA_ERROR("negative albedo");
        return 0;
    }
    return orsa::FromUnits(1329,orsa::Unit::KM)*pow(10,-0.2*H)/sqrt(p);
}

// ConstantZRotation

/* 
   bool ConstantZRotation::update(const orsa::Time & t) {
   const double _phi = _phi0 + _omega*(t-_t0).get_d();
   orsa::Matrix localMatrix = Matrix::identity();
   // same rotation as Attitude::localToGlobal(t)
   localMatrix.rotZ(_phi);
   _m = localMatrix;
   return true;
   }
*/


void orsa::principalAxis(orsa::Matrix & genericToPrincipal,
                         orsa::Matrix & principalInertiaMatrix,
                         const orsa::Matrix & inertiaMatrix) {
  
    const orsa::Matrix & I = inertiaMatrix;
  
    double data[] = { 
        I.getM11(), I.getM12(), I.getM13(), 
        I.getM21(), I.getM22(), I.getM23(), 
        I.getM31(), I.getM32(), I.getM33()};
  
    gsl_matrix_view m = gsl_matrix_view_array (data, 3, 3);
  
    gsl_vector * eval = gsl_vector_alloc (3);
  
    gsl_matrix * evec = gsl_matrix_alloc (3, 3);
  
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
  
    gsl_eigen_symmv (&m.matrix, eval, evec, w);
  
    gsl_eigen_symmv_free (w);
  
    gsl_eigen_symmv_sort (eval, evec, 
                          GSL_EIGEN_SORT_ABS_ASC);
  
    const orsa::Matrix principalToGeneric(gsl_matrix_get(evec,0,0),gsl_matrix_get(evec,0,1),gsl_matrix_get(evec,0,2),
                                          gsl_matrix_get(evec,1,0),gsl_matrix_get(evec,1,1),gsl_matrix_get(evec,1,2),
                                          gsl_matrix_get(evec,2,0),gsl_matrix_get(evec,2,1),gsl_matrix_get(evec,2,2));
  
    orsa::Matrix::transpose(principalToGeneric,genericToPrincipal);
  
    principalInertiaMatrix = genericToPrincipal*I*principalToGeneric;
  
    gsl_vector_free (eval);
    gsl_matrix_free (evec);
}

/***/

orsa::GlobalRNG * orsa::GlobalRNG::_instance = 0;

orsa::Cache<int> orsa::GlobalRNG::randomSeed;

/***/

const int RandomPointsInShape::maxRandomSeed = 1000000000;

RandomPointsInShape::RandomPointsInShape(const orsa::Shape * s,
                                         const orsa::MassDistribution * massDistribution,
                                         const size_t & samplePoints,
                                         const bool & localStoreVector) :
    osg::Referenced(),
    shape(s),
    md(massDistribution),
    size(samplePoints),
    saveVector(localStoreVector),
    randomSeed(orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform_int(maxRandomSeed)) {
    reset();
    //
    in.clear();
    density.clear();
    vec.clear();
    //
    in.resize(size);
    density.resize(size);
    if (saveVector) {
        vec.resize(size);
    }
    //
    const Box boundingBox = shape->boundingBox();
    //
    if (0) {
        //debug
        ORSA_DEBUG("box: %f x %f x %f [km]",
                   orsa::FromUnits(boundingBox.getXMax()-boundingBox.getXMin(),orsa::Unit::KM,-1),
                   orsa::FromUnits(boundingBox.getYMax()-boundingBox.getYMin(),orsa::Unit::KM,-1),
                   orsa::FromUnits(boundingBox.getZMax()-boundingBox.getZMin(),orsa::Unit::KM,-1));
    }
    //
    counter = 0;
    numInside = 0;
    while (counter < size) {
        const orsa::Vector v = __randomVectorUtil(rng.get(),boundingBox);
        in[counter] = shape->isInside(v);
        if (md.get() != 0) {
            density[counter] = md->density(v);
        }
        if (saveVector) {
            vec[counter] = v;
        }
        if (in[counter]) {
            ++numInside;
        }
        ++counter;
        
        // ORSA_DEBUG("progress: %i/%i/%i",numInside,counter,size);
    }
    reset();
}

class TriShapeInternals {
public:
    std::vector<orsa::Vector> vertexVector;
    orsa::Cache<double> volume;
    orsa::Cache<double> volumeCumulativeFraction;
};

class TriShapeInternals_LessThan {
public:
    bool operator() (const TriShapeInternals & left, const TriShapeInternals & right) {
        return left.volumeCumulativeFraction < right.volumeCumulativeFraction;
    }
    bool operator() (const TriShapeInternals & left, const double & right) {
        return left.volumeCumulativeFraction < right;
    }
    bool operator() (const double & left, const TriShapeInternals & right) {
        return left < right.volumeCumulativeFraction;
    }
};

RandomPointsInShape::RandomPointsInShape(const orsa::TriShape * s,
                                         const orsa::MassDistribution * massDistribution,
                                         const size_t & samplePoints,
                                         const bool & /* localStoreVector */) :
    osg::Referenced(),
    shape(s),
    md(massDistribution),
    size(samplePoints),
    saveVector(true),
    randomSeed(orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform_int(maxRandomSeed)) {
    reset();
    //
    in.clear();
    density.clear();
    vec.clear();
    //
    in.resize(size);
    density.resize(size);
    if (saveVector) {
        vec.resize(size);
    }
    
#warning if shape is strongly concave and a simplex covers volume outside the body shape, then the results are incorrect (including volume computations...)
    ORSA_DEBUG("warning: if shape is strongly concave and a simplex covers volume outside the body shape, then the results are incorrect (including volume computations...)");
    
    osg::ref_ptr<const orsa::TriShape> triShape = s;
    
    const orsa::TriShape::VertexVector & vv = triShape->getVertexVector();
    const orsa::TriShape::FaceVector   & fv = triShape->getFaceVector();
    std::vector<TriShapeInternals> aux;
    aux.resize(fv.size());
    for (size_t fi=0; fi<fv.size(); ++fi) {
        aux[fi].vertexVector.resize(4);
        aux[fi].vertexVector[0] = orsa::Vector(0,0,0);
        aux[fi].vertexVector[1] = vv[fv[fi].i()];
        aux[fi].vertexVector[2] = vv[fv[fi].j()];
        aux[fi].vertexVector[3] = vv[fv[fi].k()];
        aux[fi].volume =
            fabs((aux[fi].vertexVector[1]-aux[fi].vertexVector[0]) *
                 orsa::externalProduct(aux[fi].vertexVector[2]-aux[fi].vertexVector[0],
                                       aux[fi].vertexVector[3]-aux[fi].vertexVector[0]) / 6.0);
    }
    
    double totalVolume = 0.0;
    for (size_t fi=0; fi<fv.size(); ++fi) {
        totalVolume += aux[fi].volume;
    }
    
    // std::vector<double> volumeCumulativeFraction;
    // volumeCumulativeFraction.resize(fv.size());
    for (size_t fi=0; fi<fv.size(); ++fi) {
        aux[fi].volumeCumulativeFraction =  aux[fi].volume/totalVolume;
        if (fi != 0) aux[fi].volumeCumulativeFraction += aux[fi-1].volumeCumulativeFraction;
        
        // ORSA_DEBUG("volumeCumulativeFraction[%09i] = %12.9f",fi,(*aux[fi].volumeCumulativeFraction));
    }
    
    
    // set numInside so that the volume obtained from RPIS is correct
    numInside = size*totalVolume/shape->boundingBox().volume();
    
    //    Algorithm: 
    //    C. Rocchini and P. Cignoni.
    //    Generating random points in a tetrahedron.
    //    Journal of graphics tools, 5(4):9-12, 2000
    
    for (size_t counter=0; counter<size; ++counter) {
        if (counter<numInside) {
            in[counter] = true;
            
            double s = rng->gsl_rng_uniform();
            double t = rng->gsl_rng_uniform();
            double u = rng->gsl_rng_uniform();
            if (s+t>1.0) { // cut'n fold the cube into a prism
                s = 1.0 - s;
                t = 1.0 - t;
            }
            if (t+u>1.0) { // cut'n fold the prism into a tetrahedron
                double tmp = u;
                u = 1.0 - s - t;
                t = 1.0 - tmp;
            } else if (s+t+u>1.0) {
                double tmp = u;
                u = s + t + u - 1.0;
                s = 1 - t - tmp;
            }
            double a=1-s-t-u; // a,s,t,u are the barycentric coordinates of the random point.

            const double volumeCumulativeFractionSelector = rng->gsl_rng_uniform();

            std::vector<TriShapeInternals>::const_iterator it = lower_bound(aux.begin(),
                                                                            aux.end(),
                                                                            volumeCumulativeFractionSelector,
                                                                            TriShapeInternals_LessThan());
            
            // ORSA_DEBUG("val: %g lb: %g",volumeCumulativeFractionSelector,(*(*it).volumeCumulativeFraction));
            
            const orsa::Vector v =
                (*it).vertexVector[0]*a+
                (*it).vertexVector[1]*s+
                (*it).vertexVector[2]*t+
                (*it).vertexVector[3]*u;

            if (md.get() != 0) {
                density[counter] = md->density(v);
            }
            if (saveVector) {
                vec[counter] = v;
            }
            
        } else {
            in[counter] = false;
            
        }        
    }
    
    reset();
}

bool RandomPointsInShape::get(orsa::Vector & v, double & pointDensity) const {
    if (md.get() == 0) {
        ORSA_DEBUG("warning: mass distribution is not set, so the density is undefined, and this method will fail");
        ORSA_DEBUG("you should call instead get(v);");
        return false;
    }
    while (counter < size) {
        if (!saveVector) {
            // must sample, to keep the counter and the RNG in sync
            v = __randomVectorUtil(rng,shape->boundingBox());
        }
        if (in[counter]) {
            pointDensity = density[counter];
            if (saveVector) {
                v = vec[counter];
            }
            ++counter;
            return true;
        } else {
            ++counter;
        }
    }
    return false;
}

bool RandomPointsInShape::get(orsa::Vector & v) const {
    while (counter < size) {
        if (!saveVector) {
            // must sample, to keep the counter and the RNG in sync
            v = __randomVectorUtil(rng,shape->boundingBox());
        }
        if (in[counter]) {
            if (saveVector) {
                v = vec[counter];
            }
            ++counter;
            return true;
        } else {
            ++counter;
        }
    }
    return false;
}

void RandomPointsInShape::updateMassDistribution(const orsa::MassDistribution * massDistribution) {
    md = massDistribution;
    reset();
    orsa::Vector v;
    while (counter < size) {
        if (saveVector) {
            v = vec[counter];
        } else {
            v = __randomVectorUtil(rng,shape->boundingBox());
        }
        if (in[counter]) {
            density[counter] = md->density(v);
        }
        ++counter;
    }
}

orsa::Vector RandomPointsInShape::__randomVectorUtil(const orsa::RNG * rng, const Box & boundingBox) {
    return Vector(boundingBox.getXMin()+(boundingBox.getXMax()-boundingBox.getXMin())*rng->gsl_rng_uniform(),
                  boundingBox.getYMin()+(boundingBox.getYMax()-boundingBox.getYMin())*rng->gsl_rng_uniform(),
                  boundingBox.getZMin()+(boundingBox.getZMax()-boundingBox.getZMin())*rng->gsl_rng_uniform());
}

/***/

double orsa::volume(const orsa::RandomPointsInShape * randomPointsInShape) {
    return ((randomPointsInShape->shape->boundingBox().volume()*randomPointsInShape->pointsInside())/randomPointsInShape->size);
}

double orsa::mass(const orsa::RandomPointsInShape * randomPointsInShape) {
    osg::ref_ptr<orsa::Statistic<double> > stat = new orsa::Statistic<double>;
    orsa::Vector v;
    double density;
    stat->reset();
    randomPointsInShape->reset();
    while (randomPointsInShape->get(v,density)) {
        if (density > 0) {
            stat->insert(density);
        }
    }
    if (stat->entries() == 0) {
        ORSA_DEBUG("no points with positive density, returning...");
        return 0.0;
    } else {
        return stat->average()*orsa::volume(randomPointsInShape);
    }
}

orsa::Vector orsa::centerOfMass(const orsa::RandomPointsInShape * randomPointsInShape) {
    // const orsa::MassDistribution * massDistribution) {
    
    /* osg::ref_ptr<orsa::WeightedStatistic<double> > stat_CMx = new orsa::WeightedStatistic<double>;
       osg::ref_ptr<orsa::WeightedStatistic<double> > stat_CMy = new orsa::WeightedStatistic<double>;
       osg::ref_ptr<orsa::WeightedStatistic<double> > stat_CMz = new orsa::WeightedStatistic<double>;
     
       orsa::Vector v;
       randomPointsInShape->reset();
       while (randomPointsInShape->get(v)) {
       const double density = massDistribution->density(v);
       if (density > 0) {
       // stat_M->insert(density);
       //
       stat_CMx->insert(v.getX(),density);
       stat_CMy->insert(v.getY(),density);
       stat_CMz->insert(v.getZ(),density);
       }
       }
     
       const orsa::Vector center_of_mass(stat_CMx->average(),
       stat_CMy->average(),
       stat_CMz->average());
    */
  
    // one cartesian component by one, to save memory when a large number of points is used
    osg::ref_ptr<orsa::WeightedStatistic<double> > stat = new orsa::WeightedStatistic<double>;
  
    orsa::Vector v;
    double density;
    orsa::Vector center_of_mass;
    orsa::Vector center_of_mass_uncertainty;
  
    // x
    stat->reset();
    randomPointsInShape->reset();
    while (randomPointsInShape->get(v,density)) {
        // const double density = massDistribution->density(v);
        if (density > 0) {
            stat->insert(v.getX(),density);
        }
    }
    //
    if (stat->entries() == 0) {
        ORSA_DEBUG("no points with positive density, returning...");
        return orsa::Vector(0,0,0);
    }
    //
    center_of_mass.setX(stat->average());
    center_of_mass_uncertainty.setX(stat->averageError());
    // y
    stat->reset();
    randomPointsInShape->reset();
    while (randomPointsInShape->get(v,density)) {
        // const double density = massDistribution->density(v);
        if (density > 0) {
            stat->insert(v.getY(),density);
        }
    }
    center_of_mass.setY(stat->average());
    center_of_mass_uncertainty.setY(stat->averageError());
    // z
    stat->reset();
    randomPointsInShape->reset();
    while (randomPointsInShape->get(v,density)) {
        // const double density = massDistribution->density(v);
        if (density > 0) {
            stat->insert(v.getZ(),density);
        }
    }
    center_of_mass.setZ(stat->average());
    center_of_mass_uncertainty.setZ(stat->averageError());
  
    if (0) {
        // debug output
        ORSA_DEBUG("cm.x: %14.6e +/- %14.6e",
                   center_of_mass.getX(),
                   center_of_mass_uncertainty.getX());
        ORSA_DEBUG("cm.y: %14.6e +/- %14.6e",
                   center_of_mass.getY(),
                   center_of_mass_uncertainty.getY());
        ORSA_DEBUG("cm.z: %14.6e +/- %14.6e",
                   center_of_mass.getZ(),
                   center_of_mass_uncertainty.getZ());
    }
  
    return center_of_mass;
}

void orsa::diagonalizedInertiaMatrix(orsa::Matrix & shapeToLocal,
                                     orsa::Matrix & localToShape,
                                     orsa::Matrix & inertiaMatrix,
                                     const orsa::Vector & centerOfMass,
                                     const orsa::RandomPointsInShape * randomPointsInShape) {
    // const orsa::MassDistribution * massDistribution) {
    
    { 
        osg::ref_ptr<orsa::Statistic<double> > stat_Ixx = new orsa::Statistic<double>;
        osg::ref_ptr<orsa::Statistic<double> > stat_Iyy = new orsa::Statistic<double>;
        osg::ref_ptr<orsa::Statistic<double> > stat_Izz = new orsa::Statistic<double>;
        osg::ref_ptr<orsa::Statistic<double> > stat_Ixy = new orsa::Statistic<double>;
        osg::ref_ptr<orsa::Statistic<double> > stat_Ixz = new orsa::Statistic<double>;
        osg::ref_ptr<orsa::Statistic<double> > stat_Iyz = new orsa::Statistic<double>;
    
        orsa::Vector v;
        double density;
        randomPointsInShape->reset();
        while (randomPointsInShape->get(v,density)) {
            // const double density = massDistribution->density(v);
            if (density > 0) {
                v -= centerOfMass;
                // no rotation here, work in shape coordinates
                //
                stat_Ixx->insert(density*(v.getY()*v.getY()+v.getZ()*v.getZ()));
                stat_Iyy->insert(density*(v.getX()*v.getX()+v.getZ()*v.getZ()));
                stat_Izz->insert(density*(v.getX()*v.getX()+v.getY()*v.getY()));
                //
                stat_Ixy->insert(density*(-v.getX()*v.getY()));
                stat_Ixz->insert(density*(-v.getX()*v.getZ()));
                stat_Iyz->insert(density*(-v.getY()*v.getZ()));
            }
        }
        
        inertiaMatrix.set(stat_Ixx->average(),
                          stat_Ixy->average(),
                          stat_Ixz->average(),
                          stat_Ixy->average(),
                          stat_Iyy->average(),
                          stat_Iyz->average(),
                          stat_Ixz->average(),
                          stat_Iyz->average(),
                          stat_Izz->average());
        
        orsa::principalAxis(shapeToLocal,
                            inertiaMatrix,
                            inertiaMatrix);
        
        // correct, in case some of the axes got reflected
        {
            const double rot_11 = ((shapeToLocal*orsa::Vector(1,0,0))*orsa::Vector(1,0,0) < 0) ? -1 : 1;
            const double rot_22 = ((shapeToLocal*orsa::Vector(0,1,0))*orsa::Vector(0,1,0) < 0) ? -1 : 1;
            const double rot_33 = ((shapeToLocal*orsa::Vector(0,0,1))*orsa::Vector(0,0,1) < 0) ? -1 : 1;
            const orsa::Matrix rot(rot_11,0,0,
                                   0,rot_22,0,
                                   0,0,rot_33);
            shapeToLocal = rot*shapeToLocal;
        }
    
        // update localToShape as well
        orsa::Matrix::transpose(shapeToLocal,localToShape);
    
    }
  
}


orsa::PaulMoment * orsa::computePaulMoment(const unsigned int order,
                                           const orsa::Matrix & shapeToLocal,
                                           const orsa::Matrix & /* localToShape */,	
                                           const orsa::Vector & centerOfMass,
                                           const orsa::RandomPointsInShape * randomPointsInShape) {
    // const orsa::MassDistribution * massDistribution) {
    
    orsa::PaulMoment * pm = new orsa::PaulMoment(order);
  
    const unsigned int order_plus_one = order+1;
  
    osg::ref_ptr< orsa::WeightedStatistic<double> > stat = new orsa::WeightedStatistic<double>;
    
    orsa::Vector v;
    double density;
    
    for (unsigned int sum=0; sum<order_plus_one; ++sum) {
        for (unsigned int i=0; i<order_plus_one; ++i) {
            for (unsigned int j=0; j<order_plus_one; ++j) {
                for (unsigned int k=0; k<order_plus_one; ++k) {
                    if (i+j+k==sum) {
	    
                        stat->reset();
	    
                        randomPointsInShape->reset();
                        while (randomPointsInShape->get(v,density)) {
                            // const double density = massDistribution->density(v);
                            if (density > 0) {
                                v -= centerOfMass;
                                v = shapeToLocal*v;
                                stat->insert(int_pow(v.getX(),i)*
                                             int_pow(v.getY(),j)*
                                             int_pow(v.getZ(),k),
                                             density);
                            }
                        }
	    
                        // save!
                        pm->setM(stat->average(),i,j,k);
                        pm->setM_uncertainty(stat->averageError(),i,j,k);
	    
                        if (0) {
                            // debug output
                            const double  M = stat->average();
                            const double dM = stat->averageError();
                            //
                            /* ORSA_DEBUG("M[%i][%i][%i] = %g +/- %g",
                               i,j,k,M,dM);
                            */
                            //
                            const double largest = std::max(fabs(M),fabs(dM));
                            //
                            const int    p10 = floor(log10(largest));
                            const double d10 = pow(10,p10);
                            //
                            const double  M10 =  M/d10;
                            const double dM10 = dM/d10;
                            //
                            ORSA_DEBUG("M[%i][%i][%i] = (%+f +/- %f) x 10^%i",
                                       i,j,k,M10,dM10,p10);
                        }
                        
                    }
                }
            }
        }
    }
  
    return pm;
}


void orsa::bodyInertialComputations(double & volume,
                                    orsa::Vector & centerOfMass,
                                    orsa::Matrix & shapeToLocal,
                                    orsa::Matrix & localToShape,
                                    orsa::Matrix & inertiaMatrix,
                                    orsa::PaulMoment * * paulMoment,
                                    const unsigned int order,
                                    const orsa::Shape * shape,
                                    const orsa::MassDistribution * massDistribution,
                                    const unsigned int N,
                                    const bool & localStoreVector) {
    
    osg::ref_ptr<RandomPointsInShape> randomPointsInShape = new RandomPointsInShape(shape,
                                                                                    massDistribution,
                                                                                    N,
                                                                                    localStoreVector);
    
    volume = orsa::volume(randomPointsInShape);
    
    centerOfMass = orsa::centerOfMass(randomPointsInShape);
    
    orsa::diagonalizedInertiaMatrix(shapeToLocal,
                                    localToShape,
                                    inertiaMatrix,
                                    centerOfMass,
                                    randomPointsInShape);
    
    (*paulMoment) = orsa::computePaulMoment(order,
                                            shapeToLocal,
                                            localToShape,
                                            centerOfMass,
                                            randomPointsInShape);
}
