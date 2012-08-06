#include <orsa/massDistribution.h>
#include <orsa/matrix.h>
#include <orsa/paulMoment.h>
#include <orsa/util.h>
#include <orsa/vector.h>
#include <orsa/multifit.h>
#include <orsa/chebyshev.h> 

#include <orsaPDS/RadioScienceGravity.h>
#include "CubicChebyshevMassDistribution.h"
#include "CCMD2SH.h"
#include "simplex.h"
#include "SH2ijk.h"
#include "penalty.h"
#include "vesta.h"
#include "gaskell.h"

#include <qd/dd_real.h>
#include <qd/qd_real.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_siman.h>

using namespace orsa;

/*** CHOOSE ONE ***/
// typedef double simplex_T;
// typedef mpf_class simplex_T;
typedef dd_real simplex_T;
// typedef qd_real simplex_T;

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;

/*** CHOOSE ONE ***/
// typedef double T;
// typedef mpf_class T;
typedef dd_real SH_T;
// typedef qd_real T;

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SHIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SHIntegration<T>::index4Table;

//// custom mass distibutions

class ThreeComponentsMassDistribution : public MassDistribution {
public:
    ThreeComponentsMassDistribution(const orsa::Vector & coreCenter,
                                    const double & coreRx,
                                    const double & coreRy,
                                    const double & coreRz,
                                    const double & coreDensity,
                                    const double & coreMantleInterfaceThickness,
                                    const double & mantleDensity,
                                    const double & mantleCrustInterfaceThickness,
                                    const double & crustDensity,
                                    const double & crustRx,
                                    const double & crustRy,
                                    const double & crustRz) :
        MassDistribution(),
        c0(coreCenter),
        cRx(coreRx),
        cRy(coreRy),
        cRz(coreRz),
        cD(coreDensity),
        cmT(coreMantleInterfaceThickness),
        mD(mantleDensity),
        muT(mantleCrustInterfaceThickness),
        uD(crustDensity),
        uRx(crustRx),
        uRy(crustRy),
        uRz(crustRz),
        cRxp2(orsa::int_pow(cRx, 2)),
        cRyp2(orsa::int_pow(cRy, 2)),
        cRzp2(orsa::int_pow(cRz, 2)),
        uRxp2(orsa::int_pow(uRx, 2)),
        uRyp2(orsa::int_pow(uRy, 2)),
        uRzp2(orsa::int_pow(uRz, 2)),
        cRxm2(orsa::int_pow(cRx,-2)),
        cRym2(orsa::int_pow(cRy,-2)),
        cRzm2(orsa::int_pow(cRz,-2)),
        uRxm2(orsa::int_pow(uRx,-2)),
        uRym2(orsa::int_pow(uRy,-2)),
        uRzm2(orsa::int_pow(uRz,-2)) { }
protected:
    ~ThreeComponentsMassDistribution() { }
public:
    double density(const Vector & v) const {
        // all mass distribution centered in c0
        const orsa::Vector v0 = v-c0;
        const double       l0 = v0.length();
        const orsa::Vector u0 = v0.normalized();
        const double cR = sqrt(orsa::square(u0.getX())*cRxp2+
                               orsa::square(u0.getY())*cRyp2+
                               orsa::square(u0.getZ())*cRzp2);
        const double uR = sqrt(orsa::square(u0.getX())*uRxp2+
                               orsa::square(u0.getY())*uRyp2+
                               orsa::square(u0.getZ())*uRzp2);
        return (uD + (mD-uD)/(1+exp((l0-uR)/muT)) + (cD-mD)/(1+exp((l0-cR)/cmT)));
    }
public:
    const orsa::Vector c0;
    const double cRx, cRy, cRz, cD, cmT, mD, muT, uD, uRx, uRy, uRz;
    const double cRxp2, cRyp2, cRzp2, uRxp2, uRyp2, uRzp2;
    const double cRxm2, cRym2, cRzm2, uRxm2, uRym2, uRzm2;
};

// this just adds a masscon to a base massdistribution
class MassConcentrationMassDistribution : public MassDistribution {
public:
    MassConcentrationMassDistribution(const orsa::MassDistribution * baseMassDistribution,
                                      const orsa::Vector & position,
                                      const double & radius,
                                      const double & density,
                                      const double & interfaceThickness) :
        MassDistribution(),
        baseMD(baseMassDistribution),
        c0(position),
        cR(radius),
        cD(density),
        cT(interfaceThickness) { }
protected:
    ~MassConcentrationMassDistribution() { }
public:
    double density(const Vector & v) const {
        const orsa::Vector v0 = v-c0;
        const double       l0 = v0.length();
        return (baseMD->density(v) + cD/(1+exp((l0-cR)/cT)));
    }
public:
    osg::ref_ptr<const orsa::MassDistribution> baseMD;
    const orsa::Vector c0;
    const double cR, cD, cT;
};

//// ChebyshevFit3D

class ChebyshevFit3D : public orsa::Multifit {
public:
    ChebyshevFit3D(const size_t & T_degree) : 
        orsa::Multifit(),
        degree(T_degree) { }
protected:
    void singleIterationDone(const orsa::MultifitParameters *) const {
        ORSA_DEBUG("--MARK--");
    }
protected:
    const size_t degree;
    mutable std::vector< std::vector<double> >  Tx, Ty, Tz;
    void computeAllFunctionCalls(const orsa::MultifitParameters * /* par */, 
                                 const orsa::MultifitData       * data,
                                 const computeAllCallsMode        /* m */) const {
        // need to run this only once!
        static bool done=false;  
        if (!done) {
            Tx.resize(data->size());
            Ty.resize(data->size());
            Tz.resize(data->size());
            for (size_t row=0; row<data->size(); ++row) {
                orsa::ChebyshevT(Tx[row],degree,data->getD("x",row));
                orsa::ChebyshevT(Ty[row],degree,data->getD("y",row));
                orsa::ChebyshevT(Tz[row],degree,data->getD("z",row));
            }
            done=true;
        }
    }
public:
    double fun(const orsa::MultifitParameters * par, 
               const orsa::MultifitData       * /* data */,
               const unsigned int p, 
               const int          d,
               const unsigned int row) const {
        
        double f = 0.0;
        
        // const double x = data->getD("x",row);
        
        // PRECOMPUTED!
        // std::vector<double> T;
        // orsa::ChebyshevT(T,N,x);
        
        char varName[1024];
        for (size_t s=0; s<=degree; ++s) {
            for (size_t i=0; i<=degree; ++i) {
                for (size_t j=0; j<=degree-i; ++j) {
                    for (size_t k=0; k<=degree-i-j; ++k) {
                        if (i+j+k==s) {
                            sprintf(varName,"c%03i%03i%03i",i,j,k);
                            double c = par->get(varName);
                            if (p == par->index(varName)) {
                                c += d*par->getDelta(varName);
                            }
                            f += c*Tx[row][i]*Ty[row][j]*Tz[row][k];
                        }
                    }
                }
            }
        }
        
        return f;
    }
protected:
    int f_gsl(const gsl_vector * parameters, 
              void             * dataPoints, 
              gsl_vector       * f) {
    
        int retval =  orsa::Multifit::f_gsl(parameters, 
                                            dataPoints, 
                                            f);
        /* 
           const orsa::MultifitData * data = (orsa::MultifitData *) dataPoints;
           for (unsigned int j=0; j<data->size(); ++j) {
           ORSA_DEBUG("f[%02i] = %10.3f", j, gsl_vector_get(f,j));
           }
        */
        return retval;
    }
protected:
    int df_gsl (const gsl_vector * v, 
                void             * dataPoints, 
                gsl_matrix       * J) {
    
        int retval = Multifit::df_gsl(v, 
                                      dataPoints, 
                                      J);
    
        /* 
           const orsa::MultifitData * data = (orsa::MultifitData *) dataPoints;
           for (unsigned int k=0; k<_par->size(); ++k) {
           for (unsigned int j=0; j<data->size(); ++j) {
           ORSA_DEBUG("df[%02i][%02i] = %10.3f", j, k, gsl_matrix_get(J,j,k));
           }
           }
        */
        //
        return retval;
    }
  
};

/************/

// another version, using simulated annealing

/* how many points do we try before stepping */      
#define N_TRIES 100 // 100 // 200             

/* how many iterations for each T? */
#define ITERS_FIXED_T 100 // 200 // 100 // 1000 // 

/* max step size in random walk */
#define STEP_SIZE 1.0 // 1.0           

/* Boltzmann constant */
#define K 1.0                   

/* initial temperature */
// #define T_INITIAL 0.008     
#define T_INITIAL 0.010         

/* damping factor for temperature */
#warning no damping??
#define MU_T 1.000 // 1.010 // 1.003      
#define T_MIN 1.0e-5 // 2.0e-6

gsl_siman_params_t params  = {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
                              K, T_INITIAL, MU_T, T_MIN};

class SIMAN_xp {
public:
    CubicChebyshevMassDistribution::CoefficientType coeff;
    orsa::Cache<double> densityScale; // old name: bulkDensity;
    orsa::Cache<double> R0_plate;
    std::vector<orsa::Vector> rv;
    orsa::Cache<double> ref_penalty;
    std::vector<double> ref_dv;
    osg::ref_ptr<const LayerData> layerData;
};

void SIMAN_copy (void * source, void * dest) {
    SIMAN_xp * s = (SIMAN_xp *) source;
    SIMAN_xp * d = (SIMAN_xp *) dest;
    d->coeff       = s->coeff;
    d->densityScale = s->densityScale;
    d->R0_plate    = s->R0_plate;
    d->rv          = s->rv;
    d->ref_penalty = s->ref_penalty;
    d->ref_dv      = s->ref_dv;
    d->layerData   = s->layerData;
}

void * SIMAN_copy_construct (void * xp) {
    SIMAN_xp * d = new SIMAN_xp;
    SIMAN_copy(xp,d);
    return d;
}

void SIMAN_destroy (void * xp) {
    delete (SIMAN_xp *) xp;
}

double E1(void * xp) {
    
    SIMAN_xp * x = (SIMAN_xp *) xp;
    
    osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
        new CubicChebyshevMassDistribution(x->coeff,
                                           x->densityScale,
                                           x->R0_plate,
                                           x->layerData);
    
    double delta = 0.0;
    for (size_t k=0; k<x->rv.size(); ++k) {
        delta += orsa::square(massDistribution->density(x->rv[k]) - x->ref_dv[k]);
        // ORSA_DEBUG("d: %g  ref_d: %g",massDistribution->density(x->rv[k]), x->ref_dv[k]);
    }
    delta /= x->rv.size();
    delta = sqrt(delta);
    
    const double penalty =
        MassDistributionPenalty(x->rv,massDistribution.get());
    
    ORSA_DEBUG("delta: %12.6f  penalty: %12.6f   ref_penalty: %12.6f",delta,penalty,(*x->ref_penalty));
    
    if (1) {
        // another quick output...
#warning pass filename as parameter...
        CubicChebyshevMassDistributionFile::CCMDF_data data;
        data.minDensity = 0.0;
        data.maxDensity = 0.0;
        data.deltaDensity = 0.0;
        data.penalty = penalty;
        data.densityScale = x->densityScale;
        data.R0 = x->R0_plate;
        data.SH_degree = 0;
        data.coeff = x->coeff;
        data.layerData = x->layerData;
        CubicChebyshevMassDistributionFile::append(data,"MD2G.CCMDF.search.out");
    }
    
    // return (delta + 1000*orsa::square(penalty - x->ref_penalty));
    return ((1.0+penalty)*delta);
}

double M1(void * xp, void * yp) {
    SIMAN_xp * x = (SIMAN_xp *) xp;
    SIMAN_xp * y = (SIMAN_xp *) yp;
    
    double distance = 0.0;
    const size_t T_degree = x->coeff.size()-1;
    for (size_t ti=0; ti<=T_degree; ++ti) {
        for (size_t tj=0; tj<=T_degree-ti; ++tj) {
            for (size_t tk=0; tk<=T_degree-ti-tj; ++tk) {
                distance += orsa::square(y->coeff[ti][tj][tk] - x->coeff[ti][tj][tk]);
            }
        }
    }
    distance /= CubicChebyshevMassDistribution::totalSize(T_degree);
    distance = sqrt(distance);
    
    return distance;
}

void S1(const gsl_rng * r, void * xp, double step_size) {
    SIMAN_xp * x = (SIMAN_xp *) xp;
    const size_t T_degree = x->coeff.size()-1;
    for (size_t ti=0; ti<=T_degree; ++ti) {
        for (size_t tj=0; tj<=T_degree-ti; ++tj) {
            for (size_t tk=0; tk<=T_degree-ti-tj; ++tk) {
                x->coeff[ti][tj][tk] += step_size*(2*gsl_rng_uniform(r)-1)/CubicChebyshevMassDistribution::totalSize(T_degree);
            }
        }
    }    
}

void P1(void *) {
    ORSA_DEBUG("print here...");
}

/***/

int main(int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    ORSA_DEBUG("PID: %i",getpid());
    
    // QD
    unsigned int oldcw;
    fpu_fix_start(&oldcw);
    
    //ORSA_DEBUG("current mpf precision: %i",mpf_get_default_prec());
    mpf_set_default_prec(128);
    // mpf_set_default_prec(256);
    // mpf_set_default_prec(512);
    // ORSA_DEBUG("updated mpf precision: %i",mpf_get_default_prec());
    
#warning add possibility to reduce the degree...
    
#warning keep GM or change it?
    
    if ( (argc != 7) &&
         (argc != 8) ) {
        // passing CCMDF-input-file to use it as input mass distribution
        printf("Usage: %s <plate-model-file> <R0_km> <gravity-file-gravity-template-file> <output-gravity-file> <fitting-function-degree> <num-sample-points> [CCMDF-input-file]\n",argv[0]);
        exit(0);
    }
    
    const std::string plateModelFile = argv[1];
    const double plateModelR0 = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const std::string radioScienceGravityTemplateFile = argv[3];
    const std::string outputGravityFile = argv[4];
    const size_t T_degree_input = atoi(argv[5]);
    const size_t numSamplePoints = atoi(argv[6]);
    const bool have_CCMDF_file = (argc == 8);
    const std::string CCMDF_filename = (argc == 8) ? argv[7] : "";
    
    enum ALGO {
        MULTIFIT,
        ANNEALING,
        DECOMPOSITION
    };
    
    const ALGO algo = DECOMPOSITION;
    
    if (plateModelR0 <= 0.0) {
        ORSA_DEBUG("invalid input...");
        exit(0);
    }
    
    // safer over NFS
    sqlite3_vfs_register(sqlite3_vfs_find("unix-dotfile"), 1);
    
    const std::string SQLiteDBFileName = getSqliteDBFileName_simplex(plateModelFile,plateModelR0);
    
    osg::ref_ptr<GaskellPlateModel> shapeModel = new GaskellPlateModel;
    if (!shapeModel->read(plateModelFile)) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    osg::ref_ptr<SimplexIntegration<simplex_T> > si =
        new SimplexIntegration<simplex_T>(shapeModel.get(), plateModelR0, SQLiteDBFileName);
    
    osg::ref_ptr<orsaPDS::RadioScienceGravityData> gravityData = new orsaPDS::RadioScienceGravityData;
    orsaPDS::RadioScienceGravityFile::read(gravityData.get(),radioScienceGravityTemplateFile,512,1518);
    
    const double GM = gravityData->GM; 
    const double volume = si->getIntegral(0,0,0)*orsa::cube(plateModelR0);
    const double densityScale = GM/orsa::Unit::G()/volume;
    
    CubicChebyshevMassDistribution::CoefficientType densityCCC; // CCC=CubicChebyshevCoefficient
    CubicChebyshevMassDistribution::resize(densityCCC,T_degree_input);
    osg::ref_ptr<const LayerData> layerData;
    
    osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape;
    {
        const bool storeSamplePoints = true;
        randomPointsInShape =
            new orsa::RandomPointsInShape(shapeModel,
                                          0,
                                          numSamplePoints,
                                          storeSamplePoints);   
    }
    
    if (have_CCMDF_file) {
        
        CubicChebyshevMassDistributionFile::DataContainer CCMDF;
        CubicChebyshevMassDistributionFile::read(CCMDF,CCMDF_filename);
        if (CCMDF.size() == 0) {
            ORSA_DEBUG("empty CCMDF file: [%s]",CCMDF_filename.c_str());
            exit(0);
        }
        if (CCMDF.size() > 1) {
            ORSA_DEBUG("CCMDF [%s] should contain only one set of coefficients, instead of %i",CCMDF_filename.c_str(),CCMDF.size());
            /* for (size_t j=0; j<CCMDF.size(); ++j) {
               CCMDF[j].print();
               }
            */
        }
        densityCCC = CCMDF[CCMDF.size()-1].coeff;
        layerData  = CCMDF[CCMDF.size()-1].layerData;
        
    } else {
        
        // first determine the Chebyshev expansion of the mass distribution
        // #warning THIS MUST BE A PARAMETER
        const size_t T_degree = T_degree_input;
        
        // using relative density (coeff[0][0][0]=1 for constant density = bulk density)
        // CubicChebyshevMassDistribution::CoefficientType densityCCC; // CCC=CubicChebyshevCoefficient
        CubicChebyshevMassDistribution::resize(densityCCC,T_degree); 
        
        // choose mass distribution
        osg::ref_ptr<orsa::MassDistribution> massDistribution;
        //
        /* {
           const orsa::Vector coreCenter(orsa::FromUnits(0.0,orsa::Unit::KM),
           orsa::FromUnits(0.0,orsa::Unit::KM),
           orsa::FromUnits(0.0,orsa::Unit::KM));
           const double coreDensity = orsa::FromUnits(orsa::FromUnits(6.0,orsa::Unit::GRAM),orsa::Unit::CM,-3);
           // const double coreDensity = bulkDensity;
           const double mantleDensity = orsa::FromUnits(orsa::FromUnits(3.2,orsa::Unit::GRAM),orsa::Unit::CM,-3);
           // const double mantleDensity = bulkDensity;
           const double coreRadius = cbrt((3.0/(4.0*pi()))*volume*(bulkDensity-mantleDensity)/(coreDensity-mantleDensity));
           ORSA_DEBUG("coreRadius: %g [km]", orsa::FromUnits(coreRadius,orsa::Unit::KM,-1));
           massDistribution = new orsa::SphericalCorePlusMantleMassDistribution(coreCenter,
           coreRadius,
           coreDensity,
           mantleDensity);
           }
        */
        //
        {
            const double km = orsa::FromUnits(1.0,orsa::Unit::KM);
            const double gcm3 = orsa::FromUnits(orsa::FromUnits(1.0,orsa::Unit::GRAM),orsa::Unit::CM,-3);
            
            // NOTE: model densities are adjusted automatically later in order to conserve total mass
            
#warning "force" layers by using layerData?            
            // use and set layerData ??
            
            // Test
            const orsa::Vector coreCenter = orsa::Vector(0,0,0)*km;
            const double coreRx = 115*km;
            const double coreRy = 115*km;
            const double coreRz = 115*km;
            const double coreDensity = 3.4*gcm3;
            const double coreMantleInterfaceThickness = 10.0*km;
            const double mantleDensity = 3.4*gcm3;
            const double mantleCrustInterfaceThickness = 10.0*km;
            const double crustDensity = 3.4*gcm3;
            const double crustRx = 220*km;
            const double crustRy = 210*km;
            const double crustRz = 180*km;
            
            // OblateCore
            /* const orsa::Vector coreCenter = orsa::Vector(0,0,0)*km;
               const double coreRx = 130*km;
               const double coreRy = 120*km;
               const double coreRz = 100*km;
               const double coreDensity = 8.0*gcm3;
               const double coreMantleInterfaceThickness = 10.0*km;
               const double mantleDensity = 3.5*gcm3;
               const double mantleCrustInterfaceThickness = 10.0*km;
               const double crustDensity = 3.5*gcm3;
               const double crustRx = 220*km;
               const double crustRy = 210*km;
               const double crustRz = 180*km;
            */
            
            massDistribution =
                new ThreeComponentsMassDistribution(coreCenter,
                                                    coreRx,
                                                    coreRy,
                                                    coreRz,
                                                    coreDensity,
                                                    coreMantleInterfaceThickness,
                                                    mantleDensity,
                                                    mantleCrustInterfaceThickness,
                                                    crustDensity,
                                                    crustRx,
                                                    crustRy,
                                                    crustRz);
        }

        if (0) {
            // add a masscon?
            const double km = orsa::FromUnits(1.0,orsa::Unit::KM);
            const double gcm3 = orsa::FromUnits(orsa::FromUnits(1.0,orsa::Unit::GRAM),orsa::Unit::CM,-3);
            
            osg::ref_ptr<orsa::MassDistribution> baseMassDistribution = massDistribution;
            const orsa::Vector c0 = orsa::Vector(0,0,-180)*km;
            const double cR = 50*km;
            const double cD = 5.0*gcm3; // additional density, on top of baseMD density
            const double cT = 5.0*km;
            
            massDistribution =
                new MassConcentrationMassDistribution(baseMassDistribution.get(),
                                                      c0,
                                                      cR,
                                                      cD,
                                                      cT);
        }
        
        randomPointsInShape->updateMassDistribution(massDistribution.get());
        const double ref_penalty = MassDistributionPenalty(randomPointsInShape.get());
        
        ORSA_DEBUG("MD (ref) penalty: %g",ref_penalty);
        
        switch (algo) {
            case ANNEALING:
            {
                // simulated annealing 
                
                std::vector<orsa::Vector> rv;
                std::vector<double> dv;
                {
                    orsa::Vector v;
                    double density;
                    randomPointsInShape->reset();
                    while (randomPointsInShape->get(v,density)) { 
                        rv.push_back(v);
                        dv.push_back(density);
                    }
                }
                
                SIMAN_xp x0;
                x0.coeff = densityCCC;
                x0.densityScale = densityScale;
                x0.R0_plate    = plateModelR0;
                x0.rv          = rv;
                x0.ref_penalty = ref_penalty;
                x0.ref_dv      = dv;
                x0.layerData   = layerData;
                
                gsl_rng * rng = ::gsl_rng_alloc(gsl_rng_gfsr4);
                const int randomSeed = time(NULL)*getpid();
                ::gsl_rng_set(rng,randomSeed);
                ORSA_DEBUG("simulated annealing random seed: %d",randomSeed);
                
                gsl_siman_solve(rng, &x0, E1, S1, M1, P1,
                                SIMAN_copy, SIMAN_copy_construct, SIMAN_destroy,
                                0, params);
                
            }
            break;
            case MULTIFIT:
            {
                
                // multifit of mass distribution
                
                osg::ref_ptr<orsa::MultifitParameters> par = 
                    new orsa::MultifitParameters;
                char varName[1024];
                for (size_t s=0; s<=T_degree; ++s) {
                    for (size_t i=0; i<=T_degree; ++i) {
                        for (size_t j=0; j<=T_degree-i; ++j) {
                            for (size_t k=0; k<=T_degree-i-j; ++k) {
                                if (i+j+k==s) {
                                    sprintf(varName,"c%03i%03i%03i",i,j,k);
                                    par->insert(varName,0,1.0);
                                }
                            }
                        }
                    }
                }
                
                osg::ref_ptr<orsa::MultifitData> data = 
                    new orsa::MultifitData;
                data->insertVariable("x");
                data->insertVariable("y");
                data->insertVariable("z");
                orsa::Vector v;
                double density;
                const double oneOverR0 = 1.0/plateModelR0;
                randomPointsInShape->reset();
                size_t row=0;
                size_t iter=0;
                while (randomPointsInShape->get(v,density)) { 
                    data->insertD("x",row,v.getX()*oneOverR0);
                    data->insertD("y",row,v.getY()*oneOverR0);
                    data->insertD("z",row,v.getZ()*oneOverR0);
                    data->insertF(row,density/densityScale);
                    data->insertSigma(row,1.0);
                    ++row;
                    ++iter;
                }
                const size_t maxRow = --row;
                
                char logFile[1024];
                snprintf(logFile,1024,"MD2G_ChebyshevFit3D.log");
                
                osg::ref_ptr<ChebyshevFit3D> cf = new ChebyshevFit3D(T_degree);
                //
                cf->setMultifitParameters(par.get());
                cf->setMultifitData(data.get());
                //
                cf->setLogFile(logFile);
                //
                cf->run();
                
                for (size_t row=0; row<=maxRow; ++row) {
                    const double x = data->getD("x",row);
                    const double y = data->getD("y",row);
                    const double z = data->getD("z",row);
                    const double f = data->getF(row);
                    const double T = cf->fun(par.get(),
                                             data.get(),
                                             0,
                                             0,
                                             row);
                    const double err = T-f;
                    ORSA_DEBUG("FINAL: %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f",
                               x,y,z,f,T,err);
                }
                
                for (size_t s=0; s<par->totalSize(); ++s) {
                    ORSA_DEBUG("par[%03i] = [%s] = %+12.6f",s,par->name(s).c_str(),par->get(s));
                }
                
                
                for (size_t s=0; s<=T_degree; ++s) {
                    for (size_t i=0; i<=T_degree; ++i) {
                        for (size_t j=0; j<=T_degree-i; ++j) {
                            for (size_t k=0; k<=T_degree-i-j; ++k) {
                                if (i+j+k==s) {
                                    sprintf(varName,"c%03i%03i%03i",i,j,k);
                                    densityCCC[i][j][k] = par->get(varName);
                                }
                            }
                        }
                    }
                }
            }
            break;
            case DECOMPOSITION:
            {
                osg::ref_ptr<CubicChebyshevMassDistribution> CCMD =
                    CubicChebyshevMassDistributionDecomposition(massDistribution,
                                                                T_degree_input,
                                                                densityScale,
                                                                plateModelR0,
                                                                layerData);
                densityCCC = CCMD->coeff;
            }
            break;
        }
    }
    
    const size_t T_degree = densityCCC.size()-1;
    
    osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
        new CubicChebyshevMassDistribution(densityCCC,
                                           densityScale,     
                                           plateModelR0,
                                           layerData);
    randomPointsInShape->updateMassDistribution(massDistribution.get());
    
    // const double radiusCorrectionRatio = plateModelR0/gravityData->R0;

#warning not using layerData...
    
    double i1d =0.0;
    double iXd =0.0;
    double iYd =0.0;
    double iZd =0.0;
    double iXXd=0.0;
    double iXYd=0.0;
    double iXZd=0.0;
    double iYYd=0.0;
    double iYZd=0.0;
    double iZZd=0.0;
    // ti,tj,tk are the expansion of the density in terms of the cubic Chebyshev 
    for (size_t ti=0; ti<=T_degree; ++ti) {
        for (size_t tj=0; tj<=T_degree-ti; ++tj) {
            for (size_t tk=0; tk<=T_degree-ti-tj; ++tk) {
                const std::vector<mpz_class> & cTi = orsa::ChebyshevTcoeff(ti);
                const std::vector<mpz_class> & cTj = orsa::ChebyshevTcoeff(tj);
                const std::vector<mpz_class> & cTk = orsa::ChebyshevTcoeff(tk);
                // ci,cj,ck are the expansion of each Chebyshev polynomial in terms of powers of x,y,z
                for (size_t ci=0; ci<=ti; ++ci) {
                    if (cTi[ci] == 0) continue;
                    for (size_t cj=0; cj<=tj; ++cj) {
                        if (cTj[cj] == 0) continue;
                        for (size_t ck=0; ck<=tk; ++ck) {
                            if (cTk[ck] == 0) continue;
                            const double baseFactor = densityCCC[ti][tj][tk] *
                                mpz_class(cTi[ci] * cTj[cj] * cTk[ck]).get_d();
                            i1d  += baseFactor * si->getIntegral(ci,cj,ck);
                            iXd  += baseFactor * si->getIntegral(ci+1,cj,ck);
                            iYd  += baseFactor * si->getIntegral(ci,cj+1,ck);
                            iZd  += baseFactor * si->getIntegral(ci,cj,ck+1);
                            iXXd += baseFactor * si->getIntegral(ci+2,cj,ck);
                            iXYd += baseFactor * si->getIntegral(ci+1,cj+1,ck);
                            iXZd += baseFactor * si->getIntegral(ci+1,cj,ck+1);
                            iYYd += baseFactor * si->getIntegral(ci,cj+2,ck);
                            iYZd += baseFactor * si->getIntegral(ci,cj+1,ck+1);
                            iZZd += baseFactor * si->getIntegral(ci,cj,ck+2);
                        }
                    }
                }
            }
        }
    }
    const double CMx_over_plateModelR0 = iXd / i1d;
    const double CMy_over_plateModelR0 = iYd / i1d;
    const double CMz_over_plateModelR0 = iZd / i1d;
    // inertia moments, barycentric
    const double inertiaMomentXX_over_plateModelR0squared =
        ((iYYd+iZZd) - (orsa::square(CMy_over_plateModelR0)+orsa::square(CMz_over_plateModelR0))) / i1d;
    const double inertiaMomentYY_over_plateModelR0squared =
        ((iXXd+iZZd) - (orsa::square(CMx_over_plateModelR0)+orsa::square(CMz_over_plateModelR0))) / i1d;
    const double inertiaMomentZZ_over_plateModelR0squared =
        ((iXXd+iYYd) - (orsa::square(CMx_over_plateModelR0)+orsa::square(CMy_over_plateModelR0))) / i1d;
    
    ORSA_DEBUG("si->getIntegral(0,0,0): %g",si->getIntegral(0,0,0));
    ORSA_DEBUG("i1d:  %g",i1d);
    ORSA_DEBUG("iXd:  %g",iXd);
    ORSA_DEBUG("iYd:  %g",iYd);
    ORSA_DEBUG("iZd:  %g",iZd);
    ORSA_DEBUG("iXXd: %g",iXXd);
    ORSA_DEBUG("iXYd: %g",iXYd);
    ORSA_DEBUG("iXZd: %g",iXZd);
    ORSA_DEBUG("iYYd: %g",iYYd);
    ORSA_DEBUG("iYZd: %g",iYZd);
    ORSA_DEBUG("iZZd: %g",iZZd);
    ORSA_DEBUG("inertiaMomentXX_over_plateModelR0squared: %g",inertiaMomentXX_over_plateModelR0squared);
    ORSA_DEBUG("inertiaMomentYY_over_plateModelR0squared: %g",inertiaMomentYY_over_plateModelR0squared);
    ORSA_DEBUG("inertiaMomentZZ_over_plateModelR0squared: %g",inertiaMomentZZ_over_plateModelR0squared);
    
    
#warning USE THIS OR NOT??
    if (0) {
        // adjust coefficients in order to conserve total mass
        // NOTE: this changes the input densities a little (or a lot)
        //       depending on how far the generated total mass is from the nominal total mass
        
#warning not using layerData...
        
        const double correctionFactor = si->getIntegral(0,0,0)/i1d;
        for (size_t ti=0; ti<=T_degree; ++ti) {
            for (size_t tj=0; tj<=T_degree-ti; ++tj) {
                for (size_t tk=0; tk<=T_degree-ti-tj; ++tk) {
                    densityCCC[ti][tj][tk] *= correctionFactor;
                    
                    ORSA_DEBUG("coeff[%02i][%02i][%02i] = %20.12f",
                               ti,tj,tk,densityCCC[ti][tj][tk]);
                }
            }
        }
    }
    
    orsa::Cache<double> penalty;
    orsa::Cache<double> minDensity, maxDensity;
    {
        /* osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
           new CubicChebyshevMassDistribution(densityCCC,
           densityScale,    
           plateModelR0,
           layerData);
           randomPointsInShape->updateMassDistribution(massDistribution.get());
        */
        
        std::vector<orsa::Vector> rv;
        std::vector<double> dv;
        orsa::Vector v;
        double density;
        randomPointsInShape->reset();
        while (randomPointsInShape->get(v,density)) { 
            rv.push_back(v);
            dv.push_back(density);
            minDensity.setIfSmaller(density);
            maxDensity.setIfLarger(density);
        }
        penalty = MassDistributionPenalty(rv,dv,massDistribution.get());
    }
    ORSA_DEBUG("CCC penalty: %g",(*penalty));
    ORSA_DEBUG("mD: %g",(*minDensity));
    ORSA_DEBUG("MD: %g",(*maxDensity));
    
    {
        // write density,points file
        char filename[1024];
        sprintf(filename,"%s.density.points.dat",outputGravityFile.c_str());
        FILE * fp = fopen(filename,"w");
        ORSA_DEBUG("writing file [%s]",filename);
        const double gcm3 = orsa::FromUnits(orsa::FromUnits(1.0,orsa::Unit::GRAM),orsa::Unit::CM,-3);
        orsa::Vector v;
        double density;
        randomPointsInShape->reset();
        while (randomPointsInShape->get(v,density)) { 
            gmp_fprintf(fp,"%.3f\n",density/gcm3);
        }
        fclose(fp);
    }
    
    {
        // another quick output...
#warning pass filename as parameter...
        CubicChebyshevMassDistributionFile::CCMDF_data data;
#warning review all these entries
        data.minDensity = minDensity;
        data.maxDensity = maxDensity;
        data.deltaDensity = maxDensity-minDensity;
        data.penalty = penalty;
        data.densityScale = densityScale;
        data.R0 = plateModelR0;
        data.SH_degree = gravityData->degree;
        data.coeff = densityCCC;
        data.layerData = layerData;
        CubicChebyshevMassDistributionFile::write(data,"MD2G.CCMDF.out");
    }
    
    {
#warning THIS DOES NOT take into account the layerData, update!!
        
        // write barycenter file
        char filename[1024];
        sprintf(filename,"%s.inertial.dat",outputGravityFile.c_str());
        FILE * fp = fopen(filename,"w");
        ORSA_DEBUG("writing file [%s]",filename);
        const double CMx = CMx_over_plateModelR0*plateModelR0;
        const double CMy = CMy_over_plateModelR0*plateModelR0;
        const double CMz = CMz_over_plateModelR0*plateModelR0;
        gmp_fprintf(fp,"CM: %+9.3f %+9.3f %+9.3f [km]\n",
                    orsa::FromUnits(CMx,orsa::Unit::KM,-1),
                    orsa::FromUnits(CMy,orsa::Unit::KM,-1),
                    orsa::FromUnits(CMz,orsa::Unit::KM,-1));
        gmp_fprintf(fp,"Izz / (M*R0^2) = %9.6f   [R0=%g km]\n",
                    inertiaMomentZZ_over_plateModelR0squared,
                    orsa::FromUnits(plateModelR0,orsa::Unit::KM,-1));
        gmp_fprintf(fp,"Volume: %.6e [km^3]\n",orsa::FromUnits(si->getIntegral(0,0,0)*orsa::cube(plateModelR0),orsa::Unit::KM,-3));
        fclose(fp);
    }
    
#warning track precision of operations
    
    orsa::Cache<orsa::Vector> CM;
    std::vector< std::vector<mpf_class> > norm_C;
    std::vector< std::vector<mpf_class> > norm_S;
    CCMD2SH(CM,
            norm_C,
            norm_S,
            gravityData->degree,
            si.get(),
            massDistribution,
            plateModelR0,
            gravityData->R0);
    
    // only write for l>=2
    for (size_t l=2; l<=gravityData->degree; ++l) {
        for (size_t m=0; m<=l; ++m) {
            gravityData->setCoeff(orsaPDS::RadioScienceGravityData::keyC(l,m),norm_C[l][m].get_d());
            if (m!= 0) gravityData->setCoeff(orsaPDS::RadioScienceGravityData::keyS(l,m),norm_S[l][m].get_d());
        }
    }
    
    ORSA_DEBUG("writing file [%s]",outputGravityFile.c_str());
    orsaPDS::RadioScienceGravityFile::write(gravityData.get(),outputGravityFile,512,1518);
    
    return 0;
}
