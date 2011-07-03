#include "VestaInterior.h"

int main() {
    
    const size_t order = 20;
    CubicChebyshevMassDistribution::CoefficientType coeff;
    coeff.resize(order+1);
    for (size_t i=0; i<=order; ++i) {
        coeff[i].resize(order+1-i);
        for (size_t j=0; j<=order-i; ++j) {
            coeff[i][j].resize(order+1-i-j);
        }
    }
    
    // reset
    for (size_t i=0; i<=order; ++i) {
        for (size_t j=0; j<=order-i; ++j) {
            for (size_t k=0; k<=order-i-j; ++k) {
                coeff[i][j][k] = 0.0;
            }
        }
    }

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
    
    return 0;
}

