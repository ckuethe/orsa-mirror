#include "VestaInterior.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>

int main() {
    
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
    
    const size_t chebyshevDegree = 20;
    
    CubicChebyshevMassDistribution::CoefficientType coeff;
    coeff.resize(chebyshevDegree+1);
    for (size_t i=0; i<=chebyshevDegree; ++i) {
        coeff[i].resize(chebyshevDegree+1-i);
        for (size_t j=0; j<=chebyshevDegree-i; ++j) {
            coeff[i][j].resize(chebyshevDegree+1-i-j);
        }
    }
    
    // reset
    for (size_t i=0; i<=chebyshevDegree; ++i) {
        for (size_t j=0; j<=chebyshevDegree-i; ++j) {
            for (size_t k=0; k<=chebyshevDegree-i-j; ++k) {
                coeff[i][j][k] = 0.0;
            }
        }
    }
    
    const double chisq_50  = gsl_cdf_chisq_Pinv(0.50,pds->data->numberOfCoefficients);
    const double chisq_90  = gsl_cdf_chisq_Pinv(0.90,pds->data->numberOfCoefficients);
    const double chisq_95  = gsl_cdf_chisq_Pinv(0.95,pds->data->numberOfCoefficients);
    const double chisq_99  = gsl_cdf_chisq_Pinv(0.99,pds->data->numberOfCoefficients);
    //
    ORSA_DEBUG("chisq 50\%: %.2f  90\%: %.2f  95\%: %.2f  99\%: %.2f  [numberOfCoefficients=%i]",
               chisq_50,chisq_90,chisq_95,chisq_99,pds->data->numberOfCoefficients);
    
    const double g_cm3 = orsa::FromUnits(orsa::FromUnits(1,orsa::Unit::GRAM),orsa::Unit::CM,-3);
    const double maxDensity = 10.0*g_cm3;
    
    const double initialThresholdLevel = 100*chisq_99;
    const double targetThresholdLevel  = chisq_99;
    const double minAdaptiveRange      = -maxDensity;
    const double maxAdaptiveRange      =  maxDensity;
    const double intervalResidualProbability = 1.0e-10; // "1-confidence level" for this interval, different from the chisq-level
    const size_t targetSamples         = 1000;
    const size_t maxIter               = 1000000;
    
    AdaptiveIntervalVector intervalVector;
    intervalVector.resize(chebyshevDegree+1);
    for (size_t i=0; i<=chebyshevDegree; ++i) {
        intervalVector[i] = new AdaptiveIntervalType(minAdaptiveRange,
                                                     maxAdaptiveRange,
                                                     intervalResidualProbability,
                                                     initialThresholdLevel,
                                                     targetThresholdLevel,
                                                     targetSamples);
    }
    
    osg::ref_ptr<AdaptiveMonteCarloType> mc = new AdaptiveMonteCarloType;
    
    mc->run(intervalVector,maxIter);
    
    
    
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
    
    
    
    
    // free GSL stuff
    gsl_vector_free(pds_coeff);
    gsl_matrix_free(pds_covm);
    gsl_matrix_free(pds_inv_covm);
    
    return 0;
}

