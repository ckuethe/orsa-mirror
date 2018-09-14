#include "VestaSlope.h"

#include <orsa/chebyshev.h>
#include <orsa/paul.h>
#include <orsa/statistic.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>

// #include "vesta.h"
#include "gaskell.h"

#include "simplex.h"

#include "CubicChebyshevMassDistribution.h"

typedef mpfr::mpreal simplex_T;

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;

/*******/

// modified versions of RadioScienceGravityData calls, to include C10,C11,S11
unsigned int mod_gravityData_index(const orsaPDS::RadioScienceGravityData * gravityData,
                                   const QString & key) {
    unsigned int index;
    if (key == orsaPDS::RadioScienceGravityData::keyC(1,0)) {
        index = 1;
    } else if (key == orsaPDS::RadioScienceGravityData::keyC(1,1)) {
        index = 2;
    } else if (key == orsaPDS::RadioScienceGravityData::keyS(1,1)) {
        index = 3;
    } else {
        index = gravityData->index(key);
        if (index != 0) index += 3;
    }
    return index;
}

QString mod_gravityData_key(const orsaPDS::RadioScienceGravityData * gravityData,
                            const unsigned int & index) {
    if (index == 0) {
        return gravityData->key(index);
    } else if (index == 1) {
        return orsaPDS::RadioScienceGravityData::keyC(1,0);        
    } else if (index == 2) {
        return orsaPDS::RadioScienceGravityData::keyC(1,1);        
    } else if (index == 3) {
        return orsaPDS::RadioScienceGravityData::keyS(1,1);        
    } else {
        return gravityData->key(index-3);
    }
}
//
double mod_gravityData_getCoeff(const orsaPDS::RadioScienceGravityData * gravityData,
                                const QString & key) {
    double coeff;
    if ( (key == orsaPDS::RadioScienceGravityData::keyC(1,0)) ||  
         (key == orsaPDS::RadioScienceGravityData::keyC(1,1)) ||
         (key == orsaPDS::RadioScienceGravityData::keyS(1,1)) ) {
        coeff = 0.0;
    } else {
        coeff = gravityData->getCoeff(key);
    }
    return coeff;
}
//
/* double mod_gravityData_getCovar(const QString & key1, const QString & key2) {
   double covar;
   ORSA_DEBUG("complete here...");
   }
*/
//
unsigned int mod_gravityData_numberOfCoefficients(const orsaPDS::RadioScienceGravityData * gravityData) {
    return gravityData->numberOfCoefficients+3;
}

gsl_vector * mod_gravityData_getCoefficientVector(const orsaPDS::RadioScienceGravityData * gravityData) {
    gsl_vector * mu = gravityData->getCoefficientVector();
    gsl_vector * mod_mu = gsl_vector_alloc(mod_gravityData_numberOfCoefficients(gravityData));
    for (unsigned int k=0; k<mod_gravityData_numberOfCoefficients(gravityData); ++k) {
        if (k==0) {
            gsl_vector_set(mod_mu,k,gsl_vector_get(mu,k));
        } else if ( (k==1) || (k==2) || (k==3) ) {
            gsl_vector_set(mod_mu,k,0.0);
        } else {
            gsl_vector_set(mod_mu,k,gsl_vector_get(mu,k-3));
        }
    }
    return mod_mu;
}

gsl_matrix * mod_gravityData_getCovarianceMatrix(const orsaPDS::RadioScienceGravityData * gravityData) {
    gsl_matrix * covm = gravityData->getCovarianceMatrix();
    gsl_matrix * mod_covm = gsl_matrix_alloc(mod_gravityData_numberOfCoefficients(gravityData),
                                             mod_gravityData_numberOfCoefficients(gravityData));
    for (unsigned int l=0; l<mod_gravityData_numberOfCoefficients(gravityData); ++l) {
        for (unsigned int m=0; m<mod_gravityData_numberOfCoefficients(gravityData); ++m) { 
            if ((l==0) && (m==0)) {
                gsl_matrix_set(mod_covm,l,m,gsl_matrix_get(covm,l,m));
            } else if ((l==1) || (l==2) || (l==3) || (m==1) || (m==2) || (m==3)) {
                gsl_matrix_set(mod_covm,l,m,0.0);
            } else if ((l==0) && (m!=0)) {
                gsl_matrix_set(mod_covm,l,m,gsl_matrix_get(covm,l,m-3));
            } else if ((l!=0) && (m==0)) {
                gsl_matrix_set(mod_covm,l,m,gsl_matrix_get(covm,l-3,m));
            } else {
                gsl_matrix_set(mod_covm,l,m,gsl_matrix_get(covm,l-3,m-3));
            }
        }
    }
    gsl_matrix_free(covm);
    return mod_covm;   
}

gsl_matrix * mod_gravityData_getInverseCovarianceMatrix(const orsaPDS::RadioScienceGravityData * gravityData) {
    gsl_matrix * inv_covm = gravityData->getInverseCovarianceMatrix();
    gsl_matrix * mod_inv_covm = gsl_matrix_alloc(mod_gravityData_numberOfCoefficients(gravityData),
                                             mod_gravityData_numberOfCoefficients(gravityData));
    for (unsigned int l=0; l<mod_gravityData_numberOfCoefficients(gravityData); ++l) {
        for (unsigned int m=0; m<mod_gravityData_numberOfCoefficients(gravityData); ++m) { 
            if ((l==0) && (m==0)) {
                gsl_matrix_set(mod_inv_covm,l,m,gsl_matrix_get(inv_covm,l,m));
            } else if ((l==1) || (l==2) || (l==3) || (m==1) || (m==2) || (m==3)) {
                gsl_matrix_set(mod_inv_covm,l,m,0.0);
            } else if ((l==0) && (m!=0)) {
                gsl_matrix_set(mod_inv_covm,l,m,gsl_matrix_get(inv_covm,l,m-3));
            } else if ((l!=0) && (m==0)) {
                gsl_matrix_set(mod_inv_covm,l,m,gsl_matrix_get(inv_covm,l-3,m));
            } else {
                gsl_matrix_set(mod_inv_covm,l,m,gsl_matrix_get(inv_covm,l-3,m-3));
            }
        }
    }
    return mod_inv_covm;   
}

double rotationalPotential(const orsa::Vector & surfacePoint,
                           const orsa::Vector & omegaVector) {
    return (0.5*orsa::externalProduct(surfacePoint,omegaVector).lengthSquared());
}

/**********/

void gravitationalPotential_directSum(double & potential,
                                      double & min_separation,
                                      const orsa::Vector & point,
                                      const orsa::RandomPointsInShape * rpis,
                                      const double & GM,
                                      const double & smoothing_scale_squared) {
    potential = 0.0;
    min_separation = 1e20;
    if (rpis) {
        double sum_density = 0.0;
        orsa::Vector rv;
        double density;
        rpis->reset();
        while (rpis->get(rv,density)) {
            if (density <= 0.0) {
                ORSA_DEBUG("PROBLEMS... (negative or zero density) [density=%g]",density);
            }
            const orsa::Vector dr = point-rv;
            potential += density/sqrt(dr.lengthSquared()+smoothing_scale_squared);
            sum_density += density;
            if (dr.length() < min_separation) min_separation = dr.length();
        }
        potential *= GM;
        potential /= sum_density;
    }
}

void gravitationalAcceleration_directSum(orsa::Vector & acceleration,
                                         double       & min_separation,
                                         const orsa::Vector & point,
                                         const orsa::RandomPointsInShape * rpis,
                                         const double & GM,
                                         const double & smoothing_scale_squared) {
    acceleration = orsa::Vector(0,0,0);
    min_separation = 1e20;
    if (rpis) {
        double sum_density = 0.0;
        orsa::Vector rv;
        double density;
        rpis->reset();
        while (rpis->get(rv,density)) {
            const orsa::Vector dr = point-rv;
            acceleration -= dr.normalized()*density/(dr.lengthSquared()+smoothing_scale_squared);
            sum_density += density;
            if (dr.length() < min_separation) min_separation = dr.length();
        }
        acceleration *= GM;
        acceleration /= sum_density;
    }
}

//

bool geoidSearchBisection(orsa::Vector       & geoidPoint,
                          const double       & geoidPotential, /* target */
                          const orsa::Vector & surfacePlumbLine,
                          const orsa::Vector & surfacePoint,
                          const orsa::RandomPointsInShape * rpis,
                          const double       & GM,
                          const double       & smoothing_scale_squared,
                          const orsa::Vector & omegaVector,
                          const double accuracy = orsa::FromUnits(50,orsa::Unit::METER)) {
    
    /* const double rotationalPotential =
       -0.5*orsa::externalProduct(surfacePoint,omegaVector).lengthSquared();
    */
    
    const double rotPot = rotationalPotential(surfacePoint,omegaVector);
    
    // const double targetPotential = geoidPotential + rotationalPotential(surfacePoint,omegaVector);
    const double targetPotential = geoidPotential;
    
    const double accuracySquared = orsa::square(accuracy);
    
    orsa::Vector pA, pB;
    double fA, fB;
    double msA, msB;
    
    pA = surfacePoint;
    gravitationalPotential_directSum(fA,msA,pA,rpis,GM,smoothing_scale_squared);
    fA += rotPot;
    fA -= targetPotential;
    
    pB = pA;
    if (fA > 0.0) {
        do {
            pB += surfacePlumbLine.normalized()*orsa::FromUnits(10.0,orsa::Unit::KM);
            gravitationalPotential_directSum(fB,msB,pB,rpis,GM,smoothing_scale_squared);
            fB += rotPot;
            fB -= targetPotential;
        } while (fA*fB >= 0.0);
    } else {
        do {
            pB -= surfacePlumbLine.normalized()*orsa::FromUnits(10.0,orsa::Unit::KM);
            gravitationalPotential_directSum(fB,msB,pB,rpis,GM,smoothing_scale_squared);
            fB += rotPot;
            fB -= targetPotential;
        } while (fA*fB >= 0.0);
    }
    
    /* orsa::print(pA);
       orsa::print(fA);
       orsa::print(pB);
       orsa::print(fB);
    */
    
    // middle point
    orsa::Vector pM;
    double fM;
    double msM;
    
    size_t iter=0;
    bool found=false;
    do {
        ++iter;

        pM = 0.5*(pA+pB);
        gravitationalPotential_directSum(fM,msM,pM,rpis,GM,smoothing_scale_squared);
        fM += rotPot;
        fM -= targetPotential;
        
        // ORSA_DEBUG("fA: %g  fB: %g  fM: %g",fA,fB,fM);
        
        if (fM*fA <= 0.0) {
            pB = pM;
            fB = fM;
        } else {
            pA = pM;
            fA = fM;
        }
        
        if ((pB-pA).lengthSquared() < accuracySquared) {
            found=true;
            break;
        }
        
    } while (fA*fB < 0.0);
    
    /* ORSA_DEBUG("found: %i   iter: %2i   delta: %.3f [m]   sP: %.3f [km]   R: %.3f [km]   elevation: %+8.3f [km]",
       found,
       iter,
       orsa::FromUnits((pB-pA).length(),orsa::Unit::METER,-1),
       orsa::FromUnits(surfacePoint.length(),orsa::Unit::KM,-1),
       orsa::FromUnits(pM.length(),orsa::Unit::KM,-1),
       orsa::FromUnits((surfacePoint-pM)*(-surfacePlumbLine.normalized()),orsa::Unit::KM,-1));
    */
    
    if (fA*fB > 0.0) {
        ORSA_DEBUG("missed the root?");
    }
    
    geoidPoint = pM;
    
    return found;
}

orsa::Vector geoidSearchSecant_util(const double & x, const orsa::Vector & P, const orsa::Vector & u) {
    return P+u.normalized()*x;
}

bool geoidSearchSecant(orsa::Vector       & geoidPoint,
                       const double       & geoidPotential, /* target */
                       const orsa::Vector & surfacePlumbLine,
                       const orsa::Vector & surfacePoint,
                       const orsa::RandomPointsInShape * rpis,
                       const double       & GM,
                       const double       & smoothing_scale_squared,
                       const orsa::Vector & omegaVector,
                       const double accuracy = orsa::FromUnits(50,orsa::Unit::METER)) {
    
    /* const double rotationalPotential =
       -0.5*orsa::externalProduct(surfacePoint,omegaVector).lengthSquared();
    */
    
    const double rotPot = rotationalPotential(surfacePoint,omegaVector);
    
    // const double targetPotential = geoidPotential + rotationalPotential(surfacePoint,omegaVector);
    const double targetPotential = geoidPotential;
    
    ORSA_DEBUG("gP: %g   rP: %g   tP: %g",
               geoidPotential,
               rotationalPotential(surfacePoint,omegaVector),
               targetPotential);
    
    double x0, x1, x2;
    double f0, f1, f2;
    
    double ms;
    
    x0 = 0.0;
    gravitationalPotential_directSum(f0,
                                     ms,
                                     geoidSearchSecant_util(x0,surfacePoint,surfacePlumbLine),
                                     rpis,
                                     GM,
                                     smoothing_scale_squared);
    f0 += rotPot;
    f0 -= targetPotential;
    
    x1 = orsa::FromUnits(10.0,orsa::Unit::KM);
    gravitationalPotential_directSum(f1,
                                     ms,
                                     geoidSearchSecant_util(x1,surfacePoint,surfacePlumbLine),
                                     rpis,
                                     GM,
                                     smoothing_scale_squared);
    f1 += rotPot;
    f1 -= targetPotential;
    
    size_t iter=0;
    const size_t maxIter=16;
    bool converged=false;
    while(iter<maxIter) {
        
        ++iter;
        
        x2 = x1 - f1 * (x1-x0)/(f1-f0);
        
        if (fabs(x2-x1) < accuracy) {
            converged=true;          
        } else {
            gravitationalPotential_directSum(f2,
                                             ms,
                                             geoidSearchSecant_util(x2,surfacePoint,surfacePlumbLine),
                                             rpis,
                                             GM,
                                             smoothing_scale_squared);
            f2 += rotPot;
            f2 -= targetPotential;
        }
        
        if (1) {
            ORSA_DEBUG("converged: %i   iter: %2i   x: %+8.3f [km]   delta: %8.3f [km]   sP: %8.3f [km]   R: %8.3f [km]   elevation: %+8.3f [km]",
                       converged,
                       iter,
                       orsa::FromUnits(x2,orsa::Unit::KM,-1),
                       orsa::FromUnits(fabs(x2-x1),orsa::Unit::KM,-1),
                       orsa::FromUnits(surfacePoint.length(),orsa::Unit::KM,-1),
                       orsa::FromUnits(geoidSearchSecant_util(x2,surfacePoint,surfacePlumbLine).length(),orsa::Unit::KM,-1),
                       orsa::FromUnits(x2,orsa::Unit::KM,-1));
        }
        
        if (converged) {
            break;
        }
        
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f2;
    }
    
    geoidPoint = geoidSearchSecant_util(x2,surfacePoint,surfacePlumbLine);
    
    return converged;
}

//

int main(int argc, char **argv) {
    
    // orsa::GlobalRNG::randomSeed = 2117579154;
    
    orsa::Debug::instance()->initTimer();
    
    ORSA_DEBUG("PID: %i",getpid());
    
    if ( (argc != 7) &&
         (argc != 10) &&
         (argc != 11) ) {
        // printf("Usage: %s <plate-model-file> <R0_km> <RadioScienceGravityFile> <gravity-degree> [<num-sample-points> <CCMDF-input-file-1> [<CCMDF-input-file-2>]]\n",argv[0]);
        printf("Usage: %s <plate-model-file> <R0_km> <RadioScienceGravityFile> <gravity-degree> <mod-N> <mod-i> [<num-sample-points> <random-seed> <CCMDF-input-file-1> [<CCMDF-input-file-2>]]\n",argv[0]);
        exit(0);
    }   
    
    const std::string plateModelFile = argv[1];
    const double plateModelR0 = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const std::string radioScienceGravityFile = argv[3];
    const int gravityDegree = atoi(argv[4]);
    const int mod_N = atoi(argv[5]);
    const int mod_i = atoi(argv[6]);
    //
    const bool have_CCMDF_file = (argc >= 10);
    const size_t numSamplePoints = (have_CCMDF_file) ? atoi(argv[7]) : 0;
    const int randomSeed = (have_CCMDF_file) ? atoi(argv[8]) : 0;
    const std::string CCMDF_filename_1 = (argc >= 10) ? argv[9]  : "";
    const std::string CCMDF_filename_2 = (argc >= 11) ? argv[10] : "";
    
    if ( (plateModelR0 <= 0.0) ||
         (gravityDegree < 0) ||
         (mod_N < 1) ||
         (mod_i < 0) ||
         (mod_i >= mod_N) ) {
        ORSA_DEBUG("invalid input...");
        exit(0);
    }
    
    // ORSA_DEBUG("CCMDF_filename_2: [%s] == \"\" ? %i",CCMDF_filename_2.c_str(),CCMDF_filename_2=="");
    
    if (have_CCMDF_file) {
        orsa::GlobalRNG::randomSeed = randomSeed;
    }
    
    osg::ref_ptr<orsaPDS::RadioScienceGravityData> gravityData = new orsaPDS::RadioScienceGravityData;
    orsaPDS::RadioScienceGravityFile::read(gravityData.get(),radioScienceGravityFile,512,1518);
    
    osg::ref_ptr<GaskellPlateModel> shape = new GaskellPlateModel;
    if (!shape->read(plateModelFile)) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    const double volume = orsa::FromUnits(7.497e16,orsa::Unit::METER,3);
    // const double bulkDensity = gravityData->GM/orsa::Unit::G()/volume;
    
    const double km = orsa::FromUnits(1.0,orsa::Unit::KM);
    
    /* const double a = 289.0*km;
       const double b = 280.0*km;
       const double c = 229.0*km;
    */
    //
    /* const double a = 285.0*km;
       const double b = 285.0*km;
       const double c = 229.0*km;
    */
    //
    const double a = 278.2*km;
    const double b = 277.1*km;
    const double c = 231.0*km;
    
    osg::ref_ptr<orsa::EllipsoidShape> ellipsoidShape = new orsa::EllipsoidShape(a,b,c);
    
    const double T = orsa::FromUnits(5.342128,orsa::Unit::HOUR);
    const double omega = orsa::twopi()/T;
    const orsa::Vector omegaVector = orsa::Vector(0,0,omega);
    
    // computed at a radius of the volume-equivalent sphere
    // const double geoidPotentialRadius = cbrt(volume*3.0/(4*orsa::pi()));
    // geoidPotentialRadius = 285 km at the equator (including rotation effect)
    // const double geoidPotentialRadius = 285*km;
    // const double geoidPotential = gravityData->GM/geoidPotentialRadius;
    // const double geoidPotential = gravityData->GM/geoidPotentialRadius - 0.5*orsa::square(omega*geoidPotentialRadius);
    // alternative definition
    // const double geoidPotential = gravityData->GM/(230*km) - 0.333333*orsa::square(omega*(290*km));
    // const double geoidPotential = gravityData->GM/(230*km) - 0.5*orsa::square(omega*(290*km));
    // const double geoidPotential = gravityData->GM/(265*km) - 0.5*orsa::square(omega*(265*km));
    // const double geoidPotential = gravityData->GM/(255*km);
    // const double geoidPotential = gravityData->GM/(250*km);
    // const double geoidPotential = gravityData->GM/(252.4*km);
    const double geoidPotential = gravityData->GM/(251.8*km);
    // const double geoidPotential = gravityData->GM/(250.0*km);
    // ORSA_DEBUG("geoidPotentialRadius: %g   geoidPotential: %g",geoidPotentialRadius,geoidPotential);
    
    // to test direct sum for gravity
    const bool storeRandomVectors = true;
    osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape_1;
    osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape_2;
    double smoothing_scale_squared_1 = 0.0;
    double smoothing_scale_squared_2 = 0.0;
    if (have_CCMDF_file) {
        
        CubicChebyshevMassDistributionFile::DataContainer CCMDF_1;
        CubicChebyshevMassDistributionFile::read(CCMDF_1,CCMDF_filename_1);
        if (CCMDF_1.size() == 0) {
            ORSA_DEBUG("empty CCMDF file: [%s]",CCMDF_filename_1.c_str());
            exit(0);
        }
        if (CCMDF_1.size() > 1) {
            ORSA_DEBUG("CCMDF [%s] should contain only one set of coefficients.",CCMDF_filename_1.c_str());
        }
        
        const CubicChebyshevMassDistribution::CoefficientType    densityCCC_1 = CCMDF_1[CCMDF_1.size()-1].coeff;
        osg::ref_ptr<CubicChebyshevMassDistribution>  massDistribution_1 = CCMD(CCMDF_1[CCMDF_1.size()-1]);
        
        randomPointsInShape_1 =
            new orsa::RandomPointsInShape(shape,
                                          massDistribution_1.get(),
                                          numSamplePoints,
                                          storeRandomVectors);
        
        smoothing_scale_squared_1 = orsa::square(cbrt(volume/randomPointsInShape_1->pointsInside()));
        
        if (CCMDF_filename_2 != "") {
            CubicChebyshevMassDistributionFile::DataContainer CCMDF_2;
            CubicChebyshevMassDistributionFile::read(CCMDF_2,CCMDF_filename_2);

            if (CCMDF_2.size() > 0) {

                if (CCMDF_2.size() > 1) {
                    ORSA_DEBUG("CCMDF [%s] should contain only one set of coefficients.",CCMDF_filename_2.c_str());
                }
                
                const CubicChebyshevMassDistribution::CoefficientType    densityCCC_2 = CCMDF_2[CCMDF_2.size()-1].coeff;
                osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution_2 = CCMD( CCMDF_2[CCMDF_2.size()-1]);
                
#warning clone or make an independent one?
                // clone
                /* randomPointsInShape_2 = randomPointsInShape_1->clone();
                   randomPointsInShape_2->updateMassDistribution(massDistribution_2.get());
                */
                // independent
                randomPointsInShape_2 =
                    new orsa::RandomPointsInShape(shape,
                                                  massDistribution_2.get(),
                                                  numSamplePoints,
                                                  storeRandomVectors);
                
                smoothing_scale_squared_2 = orsa::square(cbrt(volume/randomPointsInShape_2->pointsInside()));
            }
        }
    }
    
    ORSA_DEBUG("smoothing_scale 1: %g [km]",sqrt(smoothing_scale_squared_1)/km);
    ORSA_DEBUG("smoothing_scale 2: %g [km]",sqrt(smoothing_scale_squared_2)/km);
    
    if (0) {
        osg::ref_ptr< orsa::Statistic<double> > stat = new orsa::Statistic<double>;
        orsa::Vector rv;
        double density;
        std::vector<orsa::Vector> vv;
        vv.reserve(numSamplePoints);
        
        vv.clear();
        randomPointsInShape_1->reset();
        while (randomPointsInShape_1->get(rv,density)) {
            vv.push_back(rv);
        }
        stat->reset();
        for (size_t j=1; j<vv.size(); ++j) {
            double closest_distance_sq = 1e99;
            for (size_t k=0; k<j; ++k) {
                if ((vv[j]-vv[k]).lengthSquared() < closest_distance_sq) {
                    closest_distance_sq = (vv[j]-vv[k]).lengthSquared();
                }
            }
            stat->insert(sqrt(closest_distance_sq));
        }
        ORSA_DEBUG("closest distance average: %g [km]   standardDeviation: %g [km]   stat size: %Zi   points inside: %i",
                   stat->average()/km,
                   stat->standardDeviation()/km,
                   stat->entries().get_mpz_t(),
                   randomPointsInShape_1->pointsInside());
        
        vv.clear();
        randomPointsInShape_2->reset();
        while (randomPointsInShape_2->get(rv,density)) {
            vv.push_back(rv);
        }
        stat->reset();
        for (size_t j=1; j<vv.size(); ++j) {
            double closest_distance_sq = 1e99;
            for (size_t k=0; k<j; ++k) {
                if ((vv[j]-vv[k]).lengthSquared() < closest_distance_sq) {
                    closest_distance_sq = (vv[j]-vv[k]).lengthSquared();
                }
            }
            stat->insert(sqrt(closest_distance_sq));
        }
        ORSA_DEBUG("closest distance average: %g [km]   standardDeviation: %g [km]   stat size: %Zi   points inside: %i",
                   stat->average()/km,
                   stat->standardDeviation()/km,
                   stat->entries().get_mpz_t(),
                   randomPointsInShape_2->pointsInside());
    }
    
    std::vector< std::vector<mpf_class> > norm_C;
    std::vector< std::vector<mpf_class> > norm_S;
    //
    norm_C.resize(gravityDegree+1);
    norm_S.resize(gravityDegree+1);
    //
    for (int l=0; l<=gravityDegree; ++l) {
        norm_C[l].resize(l+1);
        norm_S[l].resize(l+1);
        for (int m=0; m<=l; ++m) {
            if (l==0) {
                norm_C[l][m] = 1.0;
                norm_S[l][m] = 0.0;
            } else if (l==1) {
                norm_C[l][m] = 0.0;
                norm_S[l][m] = 0.0;
            } else {
                norm_C[l][m] = gravityData->getCoeff(orsaPDS::RadioScienceGravityData::keyC(l,m));
                if (m != 0) {
                    norm_S[l][m] = gravityData->getCoeff(orsaPDS::RadioScienceGravityData::keyS(l,m));
                } else {
                    norm_S[l][m] = 0.0;
                }
            }
        }
    }
    
    osg::ref_ptr<orsa::PaulMoment> paulMoment = new orsa::PaulMoment(gravityDegree);
    orsa::solve(paulMoment.get(),
                norm_C,
                norm_S,
                gravityData->R0);
    
    orsa::print(paulMoment);
    
    // double dx, dy, dz;
    orsa::Vector surfacePoint, surfaceNormal;
    const orsa::Vector origin(0,0,0);
    orsa::Vector ellipsoidPoint, ellipsoidNormal;
    
    osg::ref_ptr<orsa::PaulMoment> dummyPM = new orsa::PaulMoment(0);
    dummyPM->setM(1,0,0,0);

    char filename[1024];
    
    sprintf(filename,"slope_SH_%i_%i.xyz",mod_N,mod_i);
    FILE * fp_slope_SH_xyz = fopen(filename,"w");
    ORSA_DEBUG("writing file [%s]",filename);
    
    sprintf(filename,"slope_directSum_1_%i_%i.xyz",mod_N,mod_i);
    FILE * fp_slope_directSum_1_xyz = fopen(filename,"w");
    ORSA_DEBUG("writing file [%s]",filename);
    
    sprintf(filename,"slope_ellipsoid_%i_%i.xyz",mod_N,mod_i);
    FILE * fp_slope_ellipsoid_xyz = fopen(filename,"w");
    ORSA_DEBUG("writing file [%s]",filename);
    
    sprintf(filename,"geoid_1_%i_%i.xyz",mod_N,mod_i);
    FILE * fp_geoid_1_xyz = fopen(filename,"w");
    ORSA_DEBUG("writing file [%s]",filename);
    
    sprintf(filename,"elevation_1_%i_%i.xyz",mod_N,mod_i);
    FILE * fp_elevation_1_xyz = fopen(filename,"w");
    ORSA_DEBUG("writing file [%s]",filename);
    
    FILE * fp_slope_directSum_2_xyz = fopen(filename,"w");
    FILE * fp_geoid_2_xyz = fopen(filename,"w");
    FILE * fp_elevation_2_xyz = fopen(filename,"w");
    
    if (CCMDF_filename_2 != "") {
        
        sprintf(filename,"slope_directSum_2_%i_%i.xyz",mod_N,mod_i);
        fp_slope_directSum_2_xyz = fopen(filename,"w");
        ORSA_DEBUG("writing file [%s]",filename);
        
        sprintf(filename,"geoid_2_%i_%i.xyz",mod_N,mod_i);
        fp_geoid_2_xyz = fopen(filename,"w");
        ORSA_DEBUG("writing file [%s]",filename);
        
        sprintf(filename,"elevation_2_%i_%i.xyz",mod_N,mod_i);
        fp_elevation_2_xyz = fopen(filename,"w");
        ORSA_DEBUG("writing file [%s]",filename);
    }
    
    sprintf(filename,"slope_support_%i_%i.dat",mod_N,mod_i);
    FILE * fp_slope_support_dat = fopen(filename,"w");
    ORSA_DEBUG("writing file [%s]",filename);
    
    sprintf(filename,"arc_%i_%i.dat",mod_N,mod_i);
    FILE * fp_arc_dat = fopen(filename,"w");
    ORSA_DEBUG("writing file [%s]",filename);
    
    sprintf(filename,"arc_segment_%i_%i.dat",mod_N,mod_i);
    FILE * fp_arc_segment_dat = fopen(filename,"w");
    ORSA_DEBUG("writing file [%s]",filename);
    
    const GaskellPlateModel::VertexVector & vertexVector = shape->getVertexVector();
    size_t vertexVectorIndex = 0;
    
    // arc sampling extremes (when used...)
    // test RS
    /* const double arc_lon0 =  33.0*orsa::degToRad();
       const double arc_lat0 = -25.0*orsa::degToRad();
       const double arc_lon1 = 235.0*orsa::degToRad();
       const double arc_lat1 = -25.0*orsa::degToRad();
    */
    // good RS
    const double arc_lon0 =  37.0*orsa::degToRad();
    const double arc_lat0 = -20.0*orsa::degToRad();
    const double arc_lon1 = 230.0*orsa::degToRad();
    const double arc_lat1 = -20.0*orsa::degToRad();
    // snowman
    /* const double arc_lon0 = 180.0*orsa::degToRad();
       const double arc_lat0 =   0.0*orsa::degToRad();
       const double arc_lon1 = 215.0*orsa::degToRad();
       const double arc_lat1 =  25.0*orsa::degToRad();
    */
    // 
    double s_arc_lon, c_arc_lon;
    double s_arc_lat, c_arc_lat;
   	// sincos(arc_lon0,&s_arc_lon,&c_arc_lon);
	s_arc_lon = sin(arc_lon0);
    c_arc_lon = cos(arc_lon0);
	// sincos(arc_lat0,&s_arc_lat,&c_arc_lat);
    s_arc_lat = sin(arc_lat0);
	c_arc_lat = cos(arc_lat0);
	const orsa::Vector u_arc_0(c_arc_lon*c_arc_lat,s_arc_lon*c_arc_lat,s_arc_lat);
    // sincos(arc_lon1,&s_arc_lon,&c_arc_lon);
    s_arc_lon = sin(arc_lon1);
	c_arc_lon = cos(arc_lon1);
	// sincos(arc_lat1,&s_arc_lat,&c_arc_lat);
    s_arc_lat = sin(arc_lat1);
	c_arc_lat = cos(arc_lat1);
	const orsa::Vector u_arc_1(c_arc_lon*c_arc_lat,s_arc_lon*c_arc_lat,s_arc_lat);
    //
    const double arc_length = acos(u_arc_0*u_arc_1);
    const double sin_arc_length = sin(arc_length);
    ORSA_DEBUG("arc length: %g [deg]",arc_length*orsa::radToDeg());
    int arc_iter = 0;
    const int num_arc_segment = ceil(arc_length*orsa::radToDeg()); // approx. one per deg
    //
    const orsa::Vector u_arc_center  = (u_arc_0+u_arc_1).normalized();
    const orsa::Vector u_arc_ortho   = orsa::externalProduct(u_arc_0,u_arc_1).normalized();
    const orsa::Vector u_arc_tangent = orsa::externalProduct(u_arc_center,u_arc_ortho).normalized();
    
    while (1) {
        
        // obtain surfacePoint and surfaceNormal [two methods]
        if (1) {
            
            // CHOOSE ONE, COMMENT OUT THE OTHERS...
            
            // sample points ramdomly
            /* double dx, dy, dz;
               orsa::GlobalRNG::instance()->rng()->gsl_ran_dir_3d(&dx,&dy,&dz);
               const orsa::Vector u = orsa::Vector(dx,dy,dz).normalized();
               if (!shape->rayIntersection(surfacePoint,
               surfaceNormal,
               origin,
               u,
               false)) {
               ORSA_DEBUG("problems...");
               exit(0);
               }
            */
            
            // sample radial distance at a fixed direction
            // in this case, surfacePoint and surfaceNormal are not accurate of course 
            /* const orsa::Vector u = orsa::Vector(1,1,-1).normalized();
               surfacePoint = u * (orsa::FromUnits(0.1,orsa::Unit::KM) + orsa::FromUnits(400,orsa::Unit::KM) * orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform());
               surfaceNormal = u;
            */
            
            // orsa::print(surfacePoint);
            // orsa::print(surfaceNormal);
            
            // sample points along an arc  (using slerp for uniform sampling)
            double frac = orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
            if (arc_iter <= num_arc_segment) {
                frac = (double)arc_iter/(double)num_arc_segment;
            }
            // using slerp formula
            const orsa::Vector u = (sin(frac*arc_length)/sin_arc_length)*u_arc_0 + (sin((1.0-frac)*arc_length)/sin_arc_length)*u_arc_1;
            if (!shape->rayIntersection(surfacePoint,
                                        surfaceNormal,
                                        origin,
                                        u,
                                        false)) {
                ORSA_DEBUG("problems...");
                exit(0);
            }
            ++arc_iter;
            
        } else {
#warning USE THIS IN PRODUCTION...
            // use vertex vector in shape
            if (vertexVectorIndex<vertexVector.size()) {
                if ((vertexVectorIndex%mod_N)==(size_t)mod_i) {
                    surfacePoint = vertexVector[vertexVectorIndex];
                    surfaceNormal = shape->_getVertexNormal(vertexVectorIndex); 
                    ++vertexVectorIndex;
                } else {
                    ++vertexVectorIndex;
                    continue;
                }
            } else {
                break;
            }
        }
        
        // orsa::print(surfacePoint);
        
#warning which definition for ellipsoidPoint? closest vertex OR intersection vertex? intersection of radial or plumb line?
        //
        // ellipsoidPoint  = ellipsoidShape->closestVertex(surfacePoint);
        // ellipsoidNormal = ellipsoidShape->normalVector(ellipsoidPoint);
        //
        if (!ellipsoidShape->rayIntersection(ellipsoidPoint,
                                             ellipsoidNormal,
                                             orsa::Vector(0,0,0),
                                             surfacePoint.normalized(),
                                             false)) {
            ORSA_DEBUG("problems...");
            exit(0);
        }
        
        /*
           {
           // test: is ellipsoidPoint on the ellipsoid surface?
           const orsa::Vector & eP = ellipsoidPoint;
           ORSA_DEBUG("%g",orsa::square(eP.getX()/a)+orsa::square(eP.getY()/b)+orsa::square(eP.getZ()/c));
           }
        */
        
        // rotational
        
        // #warning keep updated all definitions of surfaceRotationalPotential in the code... (maybe write function?)
        const double surfaceRotationalPotential = rotationalPotential(surfacePoint,omegaVector);
        // -0.5*orsa::externalProduct(surfacePoint,omegaVector).lengthSquared();
        
        const orsa::Vector surfaceRotationalAcceleration =
            orsa::externalProduct(omegaVector,
                                  orsa::externalProduct(surfacePoint,omegaVector));
        
        // SH
        
        const double surfaceGravitationalPotential_SH =
            gravityData->GM * orsa::Paul::gravitationalPotential(paulMoment,
                                                                 orsa::Matrix::identity(),
                                                                 dummyPM.get(),
                                                                 orsa::Matrix::identity(),
                                                                 surfacePoint);
        
        const orsa::Vector surfaceGravitationalAcceleration_SH =
            - gravityData->GM * orsa::Paul::gravitationalForce(paulMoment,
                                                               orsa::Matrix::identity(),
                                                               dummyPM.get(),
                                                               orsa::Matrix::identity(),
                                                               surfacePoint);
        
        const orsa::Vector surfaceAcceleration_SH = 
            surfaceGravitationalAcceleration_SH +
            surfaceRotationalAcceleration;
        
        // direct sum 1
        
        orsa::Vector surfaceGravitationalAcceleration_directSum_1(0,0,0);
        double min_separation_1 = 1e20;
        gravitationalAcceleration_directSum(surfaceGravitationalAcceleration_directSum_1,
                                            min_separation_1,
                                            surfacePoint,
                                            randomPointsInShape_1.get(),
                                            gravityData->GM,
                                            smoothing_scale_squared_1);
        /* 
           if (randomPointsInShape_1.get()) {
           double sum_density = 0.0;
           orsa::Vector rv;
           double density;
           randomPointsInShape_1->reset();
           while (randomPointsInShape_1->get(rv,density)) {
           const orsa::Vector dr = surfacePoint-rv;
           surfaceGravitationalAcceleration_directSum_1 -= dr.normalized()*density/(dr.lengthSquared()+smoothing_scale_squared_1);
           sum_density += density;
           if (dr.length() < min_separation_1) min_separation_1 = dr.length();
           }
           surfaceGravitationalAcceleration_directSum_1 *= gravityData->GM;
           surfaceGravitationalAcceleration_directSum_1 /= sum_density;
           }
        */
        
        // direct sum 2
        
        orsa::Vector surfaceGravitationalAcceleration_directSum_2(0,0,0);
        double min_separation_2 = 1e20;
        if (CCMDF_filename_2 != "") {
            gravitationalAcceleration_directSum(surfaceGravitationalAcceleration_directSum_2,
                                                min_separation_2,
                                                surfacePoint,
                                                randomPointsInShape_2.get(),
                                                gravityData->GM,
                                                smoothing_scale_squared_2);
        }
        
        
        const orsa::Vector surfaceAcceleration_directSum_1 = 
            surfaceGravitationalAcceleration_directSum_1 +
            surfaceRotationalAcceleration;
        
        const orsa::Vector surfaceAcceleration_directSum_2 = 
            surfaceGravitationalAcceleration_directSum_2 +
            surfaceRotationalAcceleration;
        
        // geoid
        
        // orsa::print(surfacePoint);
        
        orsa::Vector geoidPoint_1(0,0,0);
        geoidSearchSecant(geoidPoint_1,
                          geoidPotential,
                          surfaceAcceleration_directSum_1,
                          surfacePoint, // orsa::FromUnits(300,orsa::Unit::KM)*surfacePoint.normalized(),
                          randomPointsInShape_1.get(),
                          gravityData->GM,
                          smoothing_scale_squared_1,
                          omegaVector);
        
        orsa::Vector geoidPoint_2(0,0,0);
        if (CCMDF_filename_2 != "") {
            geoidSearchSecant(geoidPoint_2,
                              geoidPotential,
                              surfaceAcceleration_directSum_2,
                              surfacePoint,
                              randomPointsInShape_2.get(),
                              gravityData->GM,
                              smoothing_scale_squared_2,
                              omegaVector);
        }
        
        //
        
        const double slope_SH = acos(-surfaceNormal*surfaceAcceleration_SH.normalized());
        
        const double slope_directSum_1 = acos(-surfaceNormal*surfaceAcceleration_directSum_1.normalized());
        
        const double slope_directSum_2 = (CCMDF_filename_2 != "") ? acos(-surfaceNormal*surfaceAcceleration_directSum_2.normalized()) : 0.0;
        
        const double slope_ellipsoid = acos(surfaceNormal*ellipsoidNormal);
        
        const double geoid_1 = geoidPoint_1.length();
        // const double geoid_1 = (geoidPoint_1-ellipsoidPoint)*(-surfaceAcceleration_directSum_1.normalized());
        
        const double geoid_2 = (CCMDF_filename_2 != "") ? geoidPoint_2.length() : 0.0;
        // const double geoid_2 = (geoidPoint_2-ellipsoidPoint)*(-surfaceAcceleration_directSum_2.normalized());
        
        const double elevation_1 = (surfacePoint-geoidPoint_1)*(-surfaceAcceleration_directSum_1.normalized());
        
        const double elevation_2 = (CCMDF_filename_2 != "") ? (surfacePoint-geoidPoint_2)*(-surfaceAcceleration_directSum_2.normalized()) : 0.0;

        {
            // if (slope_SH < 89.9*orsa::degToRad()) {
            // avoid spurious points at 90.0 deg
            
            // ORSA_DEBUG("%g",orsa::radToDeg()*acos(-surfaceNormal*surfaceAcceleration.normalized()));
            
            const orsa::Vector & v = surfacePoint;
            const double lat = asin(v.getZ()/v.length());
            const double lon = fmod(orsa::twopi()+atan2(v.getY(),v.getX()),orsa::twopi());
            
            fprintf(fp_slope_SH_xyz,
                    "%g %g %g\n",
                    lon*orsa::radToDeg(),
                    lat*orsa::radToDeg(),
                    orsa::radToDeg()*slope_SH);
            fflush(fp_slope_SH_xyz);

            fprintf(fp_slope_directSum_1_xyz,
                    "%g %g %g\n",
                    lon*orsa::radToDeg(),
                    lat*orsa::radToDeg(),
                    orsa::radToDeg()*slope_directSum_1);
            fflush(fp_slope_directSum_1_xyz);
            
            fprintf(fp_slope_ellipsoid_xyz,
                    "%g %g %g\n",
                    lon*orsa::radToDeg(),
                    lat*orsa::radToDeg(),
                    orsa::radToDeg()*slope_ellipsoid);
            fflush(fp_slope_ellipsoid_xyz);
            
            fprintf(fp_geoid_1_xyz,
                    "%g %g %g\n",
                    lon*orsa::radToDeg(),
                    lat*orsa::radToDeg(),
                    geoid_1/km);
            fflush(fp_geoid_1_xyz);
            
            fprintf(fp_elevation_1_xyz,
                    "%g %g %g\n",
                    lon*orsa::radToDeg(),
                    lat*orsa::radToDeg(),
                    elevation_1/km);
            fflush(fp_elevation_1_xyz);
            
            if (CCMDF_filename_2 != "") {
                
                fprintf(fp_slope_directSum_2_xyz,
                        "%g %g %g\n",
                        lon*orsa::radToDeg(),
                        lat*orsa::radToDeg(),
                        orsa::radToDeg()*slope_directSum_2);
                fflush(fp_slope_directSum_2_xyz);
                
                fprintf(fp_geoid_2_xyz,
                        "%g %g %g\n",
                        lon*orsa::radToDeg(),
                        lat*orsa::radToDeg(),
                        geoid_2/km);
                fflush(fp_geoid_2_xyz);
                
                fprintf(fp_elevation_2_xyz,
                        "%g %g %g\n",
                        lon*orsa::radToDeg(),
                        lat*orsa::radToDeg(),
                        elevation_2/km);
                fflush(fp_elevation_2_xyz);
                
            }
            
#warning ADD geoidPoint info and elevation!
#warning if changing this format, update the \"craterslope\" program as well, which reads this file
            fprintf(fp_slope_support_dat,
                    "%g %g %g %g   %g %g %g %g   %g %g   %g %g %g   %g %g %g %g %g %g %g %g   %g %g   %g %g %g   %g %g %g   %g %g %g   %g %g %g\n",
                    surfacePoint.getX()/km,
                    surfacePoint.getY()/km,
                    surfacePoint.getZ()/km,
                    surfacePoint.length()/km,
                    //
                    ellipsoidPoint.getX()/km,
                    ellipsoidPoint.getY()/km,
                    ellipsoidPoint.getZ()/km,
                    ellipsoidPoint.length()/km,
                    //
                    lon*orsa::radToDeg(),
                    lat*orsa::radToDeg(),                    
                    //
                    orsa::radToDeg()*acos(surfaceNormal*ellipsoidNormal),
                    orsa::radToDeg()*acos(surfaceNormal*surfacePoint.normalized()),
                    orsa::radToDeg()*acos(ellipsoidNormal*ellipsoidPoint.normalized()),
                    //
                    
                    //
                    surfaceGravitationalAcceleration_SH.length(),
                    surfaceRotationalAcceleration.length(),
                    surfaceAcceleration_SH.length(),
                    gravityData->GM/surfacePoint.lengthSquared(),
                    surfaceGravitationalAcceleration_directSum_1.length(),
                    (CCMDF_filename_2 != "") ? surfaceGravitationalAcceleration_directSum_2.length() : 0.0,
                    surfaceAcceleration_directSum_1.length(),
                    (CCMDF_filename_2 != "") ? surfaceAcceleration_directSum_2.length() : 0.0,
                    //
                    min_separation_1/km,
                    min_separation_2/km,
                    //
                    orsa::radToDeg()*acos(surfaceGravitationalAcceleration_SH.normalized()*surfaceGravitationalAcceleration_directSum_1.normalized()),
                    (CCMDF_filename_2 != "") ? orsa::radToDeg()*acos(surfaceGravitationalAcceleration_SH.normalized()*surfaceGravitationalAcceleration_directSum_2.normalized()) : 0.0,
                    (CCMDF_filename_2 != "") ? orsa::radToDeg()*acos(surfaceGravitationalAcceleration_directSum_1.normalized()*surfaceGravitationalAcceleration_directSum_2.normalized()) : 0.0,
                    //
                    orsa::radToDeg()*acos(surfaceNormal*surfaceAcceleration_SH.normalized()),
                    orsa::radToDeg()*acos(surfaceNormal*surfaceAcceleration_directSum_1.normalized()),
                    (CCMDF_filename_2 != "") ? orsa::radToDeg()*acos(surfaceNormal*surfaceAcceleration_directSum_2.normalized()) : 0.0,
                    //
                    orsa::radToDeg()*acos(surfacePoint.normalized()*surfaceAcceleration_SH.normalized()),
                    orsa::radToDeg()*acos(surfacePoint.normalized()*surfaceAcceleration_directSum_1.normalized()),
                    (CCMDF_filename_2 != "") ? orsa::radToDeg()*acos(surfacePoint.normalized()*surfaceAcceleration_directSum_2.normalized()) : 0.0,
                    //
                    orsa::radToDeg()*acos(ellipsoidNormal*surfaceAcceleration_SH.normalized()),
                    orsa::radToDeg()*acos(ellipsoidNormal*surfaceAcceleration_directSum_1.normalized()),
                    (CCMDF_filename_2 != "") ? orsa::radToDeg()*acos(ellipsoidNormal*surfaceAcceleration_directSum_2.normalized()) : 0.0);
            fflush(fp_slope_support_dat);
            
            {
                // arc
                const double arc_point = acos(surfacePoint.normalized()*u_arc_0);
                const double arg = 0.5*arc_length-arc_point;
                // double s_arg, c_arg;
                // sincos(arg,&s_arg,&c_arg);
                fprintf(fp_arc_dat,
                        "%g %g   %g %g   %g %g %g   %g %g %g %g %g %g\n",
                        arc_point*orsa::radToDeg(),
                        arg*orsa::radToDeg(),
                        //
                        lon*orsa::radToDeg(),
                        lat*orsa::radToDeg(),
                        //
                        ellipsoidPoint.length()/km,
                        geoidPoint_1.length()/km,
                        surfacePoint.length()/km,
                        //
                        ellipsoidPoint*u_arc_tangent/km,
                        ellipsoidPoint*u_arc_center/km,
                        geoidPoint_1*u_arc_tangent/km,
                        geoidPoint_1*u_arc_center/km,
                        surfacePoint*u_arc_tangent/km,
                        surfacePoint*u_arc_center/km);
                fflush(fp_arc_dat);

                if (arc_iter <= num_arc_segment) {
                    // two lines for each segment entry
                    double segmentLength = orsa::FromUnits(30,orsa::Unit::KM);
                    //
                    // segmentLength *= surfaceAcceleration_directSum_1.length()/0.23; // proportional to surface acceleration
                    //
                    /* fprintf(fp_arc_segment_dat,"%g %g\n",
                       surfacePoint*u_arc_tangent/km,
                       surfacePoint*u_arc_center/km);
                       fprintf(fp_arc_segment_dat,"%g %g\n",
                       (surfacePoint+surfaceAcceleration_directSum_1.normalized()*segmentLength)*u_arc_tangent/km,
                       (surfacePoint+surfaceAcceleration_directSum_1.normalized()*segmentLength)*u_arc_center/km);
                    */
                    //
                    fprintf(fp_arc_segment_dat,"%g %g\n",
                            surfacePoint*u_arc_tangent/km,
                            surfacePoint*u_arc_center/km);
                    fprintf(fp_arc_segment_dat,"%g %g\n",
                            (geoidPoint_1+surfaceAcceleration_directSum_1.normalized()*segmentLength)*u_arc_tangent/km,
                            (geoidPoint_1+surfaceAcceleration_directSum_1.normalized()*segmentLength)*u_arc_center/km);
                    //
                    /* 
                       fprintf(fp_arc_segment_dat,"%g %g\n",
                       (geoidPoint_1-surfaceAcceleration_directSum_1.normalized()*0.5*segmentLength)*u_arc_tangent/km,
                       (geoidPoint_1-surfaceAcceleration_directSum_1.normalized()*0.5*segmentLength)*u_arc_center/km);
                       fprintf(fp_arc_segment_dat,"%g %g\n",
                       (geoidPoint_1+surfaceAcceleration_directSum_1.normalized()*segmentLength)*u_arc_tangent/km,
                       (geoidPoint_1+surfaceAcceleration_directSum_1.normalized()*segmentLength)*u_arc_center/km);
                    */
                    //
                    fflush(fp_arc_segment_dat);
                }
                
            }
            
        }
        
    }
    
    fclose(fp_slope_SH_xyz);
    fclose(fp_slope_directSum_1_xyz);
    fclose(fp_slope_ellipsoid_xyz);
    fclose(fp_geoid_1_xyz);
    fclose(fp_elevation_1_xyz);
    if (CCMDF_filename_2 != "") {
        fclose(fp_slope_directSum_2_xyz);
        fclose(fp_geoid_2_xyz);
        fclose(fp_elevation_2_xyz);
    }    
    fclose(fp_slope_support_dat);
    fclose(fp_arc_dat);
    fclose(fp_arc_segment_dat);
    
    return 0;
}

