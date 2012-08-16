#include "CubicChebyshevMassDistribution.h"

int main(int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if (argc != 4) {
        printf("Usage: %s <CCMD-densityScale> <R0_km> <CCMD-degree>\n",argv[0]);
        exit(0);
    }
    
    const double densityScale = orsa::FromUnits(orsa::FromUnits(atof(argv[1]),orsa::Unit::GRAM),orsa::Unit::CM,-3);
    const double plateModelR0 = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const size_t T_degree_input = atoi(argv[3]);
    
    CubicChebyshevMassDistribution::CoefficientType densityCCC;
    CubicChebyshevMassDistribution::resize(densityCCC,T_degree_input);

    CubicChebyshevMassDistributionFile::CCMDF_data data;
    data.minDensity   = 0.0;
    data.maxDensity   = 0.0;
    data.deltaDensity = 0.0;
    data.penalty      = 0.0;
    data.densityScale = densityScale;
    data.R0           = plateModelR0;
    data.SH_degree    = 0;
    data.coeff        = densityCCC;
    data.layerData    = 0;
    CubicChebyshevMassDistributionFile::write(data,"CCMDF.dummy.out");
    
    return 0;
}
