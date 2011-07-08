#include "VestaInteriorMonteCarlo.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>

#include "vesta.h"

int main() {
    
    // test specific cases, for debug purposes only!
    // orsa::GlobalRNG::randomSeed = -800402816;
    
    orsa::Debug::instance()->initTimer();
    
    osg::ref_ptr<orsaPDS::RadioScienceGravityFile> pds =
        new orsaPDS::RadioScienceGravityFile("JGDAWN20SIMA.DAT",512,1518);
    
    gsl_vector * pds_coeff    = pds->data->getCoefficientVector();
    gsl_matrix * pds_covm     = pds->data->getCovarianceMatrix();
    gsl_matrix * pds_inv_covm = pds->data->getInverseCovarianceMatrix();
    
    if (0) {
        // TEST ONLY
        gsl_matrix * identity = gsl_matrix_alloc(pds->data->numberOfCoefficients,pds->data->numberOfCoefficients);
        // gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,pds_covm,pds_inv_covm,0.0,identity);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,pds_inv_covm,pds_covm,0.0,identity);
        for (size_t i=0; i<pds->data->numberOfCoefficients; ++i) {
            for (size_t j=0; j<pds->data->numberOfCoefficients; ++j) {
                ORSA_DEBUG("identity[%03i][%03i] = %+20.9f",i,j,gsl_matrix_get(identity,i,j));
            }
        }
    }
    
    osg::ref_ptr<VestaShape> shape = new VestaShape;
    if (!shape->read("vesta_thomas.dat")) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    const double g_cm3 = orsa::FromUnits(orsa::FromUnits(1,orsa::Unit::GRAM),orsa::Unit::CM,-3);
    const double maxDensity = 10.0*g_cm3;
    
    const size_t chebyshevDegree = 2;
    osg::ref_ptr<AuxiliaryData> auxiliaryData = new AuxiliaryData;
    auxiliaryData->shape = shape;
    auxiliaryData->sphericalHarmonicDegree = 4; // pds->data->degree;
    auxiliaryData->chebyshevDegree = chebyshevDegree;
    auxiliaryData->R0 = pds->data->R0;
    auxiliaryData->numSamplePoints = 1000000; // MonteCarlo to determine spherical harmonics coefficients
    auxiliaryData->storeSamplePoints = true;
    auxiliaryData->intervalVectorSize = CubicChebyshevMassDistribution::totalSize(chebyshevDegree);
    auxiliaryData->pds_data = pds->data.get();
    auxiliaryData->pds_numberOfCoefficients = pds->data->numberOfCoefficients;
    auxiliaryData->pds_coeff = pds_coeff;
    auxiliaryData->pds_inv_covm = pds_inv_covm;
    
    ORSA_DEBUG("intervalVectorSize: %i",(*auxiliaryData->intervalVectorSize));

#warning UPDATE DOF if you include more vars, i,e, center of mass position, body volume...
#warning double-check this definition...
    // DOF = (sphericalHarmonicDegree+1)^2 + 1 , this last +1 is due to the GM term
    size_t chisq_DOF = (auxiliaryData->sphericalHarmonicDegree+1)*(auxiliaryData->sphericalHarmonicDegree+1)+1;
    
    const double chisq_50  = gsl_cdf_chisq_Pinv(0.50,chisq_DOF);
    const double chisq_90  = gsl_cdf_chisq_Pinv(0.90,chisq_DOF);
    const double chisq_95  = gsl_cdf_chisq_Pinv(0.95,chisq_DOF);
    const double chisq_99  = gsl_cdf_chisq_Pinv(0.99,chisq_DOF);
    //
    ORSA_DEBUG("chisq 50\%: %.2f  90\%: %.2f  95\%: %.2f  99\%: %.2f  [DOF=%i]",
               chisq_50,chisq_90,chisq_95,chisq_99,chisq_DOF);
    
    const double initialThresholdLevel = 1e20; // 100*chisq_99;
    const double targetThresholdLevel  = chisq_99;
    // const double minAdaptiveRange      = -maxDensity;
    // const double maxAdaptiveRange      =  maxDensity;
    const double intervalResidualProbability = 1.0e-10; // "1-confidence level" for this interval, different from the chisq-level
    const size_t targetSamples         = 1000;
    const size_t maxIter               = 1000000;
    
    AdaptiveIntervalVector intervalVector;
    intervalVector.resize(auxiliaryData->intervalVectorSize);
    for (unsigned int i=0; i<=auxiliaryData->chebyshevDegree; ++i) {
        for (unsigned int j=0; j<=auxiliaryData->chebyshevDegree; ++j) {
            for (unsigned int k=0; k<=auxiliaryData->chebyshevDegree; ++k) {
                if (i+j+k<=auxiliaryData->chebyshevDegree) {
                    const size_t index = CubicChebyshevMassDistribution::index(i,j,k);
                    const size_t degree = i+j+k;
#warning keep an eye on these limits
                    double minAdaptiveRange      = -maxDensity/(degree+1);
                    double maxAdaptiveRange      =  maxDensity/(degree+1);
                    if (degree==0) {
                        minAdaptiveRange = 0.0;
                    }
                    if ((i%2==1) || (j%2==1) || (k%2==1)) {
                        // odd degree, or even degree with one odd index
                        minAdaptiveRange *= 0.1;
                        maxAdaptiveRange *= 0.1;
                    }
                    intervalVector[index] =
                        new AdaptiveIntervalType(minAdaptiveRange,
                                                 maxAdaptiveRange,
                                                 intervalResidualProbability,
                                                 initialThresholdLevel,
                                                 targetThresholdLevel,
                                                 targetSamples,
                                                 auxiliaryData.get());
                }
            }
        }
    }
    
    osg::ref_ptr<AdaptiveMonteCarloType> mc = new AdaptiveMonteCarloType(auxiliaryData.get());
    
    mc->run(intervalVector,maxIter);
    
    
    /* 
       if (0) { 
       // TEST ONLY
       
       // set values
       coeff[0][0][0] = 1.0;
       coeff[0][0][1] = 0.1;
       coeff[0][1][1] = 0.2;
       coeff[2][1][3] = 0.2;
       
       osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
       new CubicChebyshevMassDistribution(coeff);
       
       {
       double x,y,z;
       for (size_t s=0; s<1000000; ++s) {
       orsa::GlobalRNG::instance()->rng()->gsl_ran_dir_3d(&x,&y,&z);
       const orsa::Vector p(x,y,z);
       const double density = massDistribution->density(p);
       ORSA_DEBUG("SAMPLE %+9.6f %+9.6f %+9.6f %+9.6f",x,y,z,density);
       }
       }
       }
    */
    
    
    // free GSL stuff
    gsl_vector_free(pds_coeff);
    gsl_matrix_free(pds_covm);
    gsl_matrix_free(pds_inv_covm);
    
    return 0;
}

