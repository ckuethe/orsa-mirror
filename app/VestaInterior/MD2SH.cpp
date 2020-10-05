#include <orsa/massDistribution.h>
#include <orsa/matrix.h>
#include <orsa/paulMoment.h>
#include <orsa/util.h>
#include <orsa/vector.h>
#include <orsa/multifit.h>
#include <orsa/chebyshev.h> 
#include <orsa/statistic.h>

#include <orsaPDS/RadioScienceGravity.h>
#include "CubicChebyshevMassDistribution.h"
#include "CCMD2SH.h"
#include "CCMD2ijk.h"
#include "simplex.h"
#include "SH2ijk.h"
#include "penalty.h"
#include "shape.h"
#include "translate_ijk.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_sf_legendre.h>

using namespace orsa;

#include "mpreal.h"

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;

/***/

typedef mpfr::mpreal F;

// calculate gravity SH expansion the "hard" way, from MonteCarlo integrals... for testing purposes mainly

int main(int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    ORSA_DEBUG("PID: %i",getpid());    
    
    if (argc != 6) {
        printf("Usage: %s <plate-model-file> <R0_km> <gravity-R0_km> <max-degree-gravity> <CCMDF-input-file>\n",argv[0]);
        exit(0);
    }
    
    const double km = orsa::FromUnits(1.0,orsa::Unit::KM);
    
    const std::string plateModelFile = argv[1];
    const double plateModelR0 = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const double gravityR0 = orsa::FromUnits(atof(argv[3]),orsa::Unit::KM);
    const int lmax = atoi(argv[4]);
    const std::string CCMDF_filename = argv[5];
    
    if (plateModelR0 <= 0.0) {
        ORSA_DEBUG("invalid input...");
        exit(0);
    }
    
    osg::ref_ptr<InputShape> shapeModel = new InputShape;
       if (!shapeModel->read(plateModelFile)) {
       ORSA_ERROR("problems encountered while reading shape file...");
       exit(0);
    }
    
    const std::string SQLiteDBFileName = getSqliteDBFileName_simplex(plateModelFile,plateModelR0);
    osg::ref_ptr<SimplexIntegration<F> > si = new SimplexIntegration<F>(shapeModel.get(), plateModelR0, SQLiteDBFileName);
    
    // osg::ref_ptr<orsaPDS::RadioScienceGravityData> gravityData = new orsaPDS::RadioScienceGravityData;
    // orsaPDS::RadioScienceGravityFile::read(gravityData.get(),radioScienceGravityTemplateFile,512,1518);
    
    // const double GM = gravityData->GM; 
    const double volume = si->getIntegral(0,0,0)*orsa::cube(plateModelR0);
    // const double densityScale = GM/orsa::Unit::G()/volume;
    // double densityScale = 1000.0;
    
    // CubicChebyshevMassDistribution::CoefficientType densityCCC; // CCC=CubicChebyshevCoefficient
    // CubicChebyshevMassDistribution::resize(densityCCC,T_degree_input);
    // osg::ref_ptr<const LayerData> layerData;
    
    /*
    osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape;
    {
        randomPointsInShape =
            new orsa::RandomPointsInShape(shapeModel,
                                          0,
                                          numSamplePoints,
                                          storeSamplePoints);   
    }
    */
    
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
    CubicChebyshevMassDistribution::CoefficientType densityCCC = CCMDF[CCMDF.size()-1].coeff;
    osg::ref_ptr<const LayerData> layerData = CCMDF[CCMDF.size()-1].layerData;
    // densityScale = CCMDF[CCMDF.size()-1].densityScale;
    
    osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
        new CubicChebyshevMassDistribution(densityCCC,
                                           // densityScale,     
                                           plateModelR0,
                                           layerData);
                                           
    /*
    std::vector< std::vector<mpf_class> > norm_C;
    std::vector< std::vector<mpf_class> > norm_S;
    norm_C.resize(lmax+1);
    norm_S.resize(lmax+1);
    for (size_t l=0; l<=lmax; ++l) {
        norm_C[l].resize(l+1);
        norm_S[l].resize(l+1);
        for (size_t m=0; m<=l; ++m) {
            norm_C[l][m]=0.0;
            norm_S[l][m]=0.0;
        }
    }
    */
    
    std::vector< std::vector< osg::ref_ptr< orsa::Statistic<double> > > > norm_C;
    std::vector< std::vector< osg::ref_ptr< orsa::Statistic<double> > > > norm_S;
    norm_C.resize(lmax+1);
    norm_S.resize(lmax+1);
    for (size_t l=0; l<=lmax; ++l) {
        norm_C[l].resize(l+1);
        norm_S[l].resize(l+1);
        for (size_t m=0; m<=l; ++m) {
            norm_C[l][m] = new orsa::Statistic<double>;
            norm_S[l][m] = new orsa::Statistic<double>;
        }
    }
    
    const bool storeSamplePoints = true;
    const int numSamplePoints=100000; // this is the incremental number, it's in a loop...
    orsa::Vector v;
    double density;
    while (1) {
        osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape = new orsa::RandomPointsInShape(shapeModel,massDistribution,numSamplePoints,storeSamplePoints);   
        while (randomPointsInShape->get(v,density)) {
            
            const orsa::Vector u = v.normalized();
            const double c_theta = u.getZ();
            const double phi = atan2(u.getY(),u.getX());
            
            double * Plm_array = new double[gsl_sf_legendre_array_n(lmax)];
            gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,lmax,c_theta,+1,Plm_array);
            
            for (size_t l=0; l<=lmax; ++l) {
                for (size_t m=0; m<=l; ++m) {
                    const double Plm = Plm_array[gsl_sf_legendre_array_index(l,m)];
                              norm_C[l][m]->insert(density*Plm*pow(v.length()/gravityR0,l)*cos(m*phi));
                    if (m!=0) norm_S[l][m]->insert(density*Plm*pow(v.length()/gravityR0,l)*sin(m*phi));
                }
            }
            delete[] Plm_array;
        }
        
        // output
        for (size_t l=0; l<=lmax; ++l) {
            for (size_t m=0; m<=l; ++m) {
                const double nf = orsa::normalization_integralToNormalizedSphericalHarmonics(l,m).get_d()/norm_C[0][0]->average(); // normalization factor
                          ORSA_DEBUG("norm_C[%02i][%02i] = %+12.9f +/- %+12.9f",l,m,nf*norm_C[l][m]->average(),nf*norm_C[l][m]->averageError());
                if (m!=0) ORSA_DEBUG("norm_S[%02i][%02i] = %+12.9f +/- %+12.9f",l,m,nf*norm_S[l][m]->average(),nf*norm_S[l][m]->averageError());
            }
        }
    }

    return 0;
}
