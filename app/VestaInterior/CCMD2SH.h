#ifndef CCMD2SH_H
#define CCMD2SH_H

#include "CubicChebyshevMassDistribution.h"
#include "simplex.h"

template <typename T>
void CCMD2SH(orsa::Vector & CM, /* note: CM is just an output variable */
             std::vector< std::vector<mpf_class> > & norm_C,
             std::vector< std::vector<mpf_class> > & norm_S,
             const size_t                          & degree,
             const SimplexIntegration<T>           * si,
             const CubicChebyshevMassDistribution::CoefficientType & coeff,
             const double                          & plateModelR0,
             const double                          & gravityDataR0);

#endif // CCMD2SH_H
