#ifndef CCMD2SH_H
#define CCMD2SH_H

#include "CubicChebyshevMassDistribution.h"
#include "simplex.h"
#include "SH2ijk.h"
#include "global_SH_epsrel.h"
#include "CCMD2ijk.h"

// note on CM: it is used as input if set, or it is computed using CCMD if unset
template <typename T>
void CCMD2SH(orsa::Cache<orsa::Vector>             & CM,
             std::vector< std::vector<mpf_class> > & norm_C,
             std::vector< std::vector<mpf_class> > & norm_S,
             // std::vector< std::vector< std::vector<double> > > & global_N,
             // std::vector< std::vector< std::vector<double> > > & translated_global_N,
             const size_t                          & SH_degree,
             const SimplexIntegration<T>           * si,
             const CubicChebyshevMassDistribution  * CCMD,
             const double                          & plateModelR0,
             const double                          & gravityDataR0) {
    
    // std::vector< std::vector< std::vector<double> > > & N = global_N;
    std::vector< std::vector< std::vector<double> > > N;
    CCMD2ijk(N,
             SH_degree,
             si,
             CCMD,
             plateModelR0);
    
    const orsa::Vector CM_over_plateModelR0 =
        (CM.isSet()) ? (CM/plateModelR0) : orsa::Vector(N[1][0][0]/N[0][0][0],
                                                        N[0][1][0]/N[0][0][0],
                                                        N[0][0][1]/N[0][0][0]);
    if (!CM.isSet()) {
        CM = CM_over_plateModelR0*plateModelR0;
    }
    
    // std::vector< std::vector< std::vector<double> > > & translated_N = translated_global_N;
    std::vector< std::vector< std::vector<double> > > translated_N;
    translate(translated_N,N,-CM_over_plateModelR0);
    if (1) {
        for (size_t degree=0; degree<=SH_degree; ++degree) {
            for (size_t ni=0; ni<=degree; ++ni) {
                for (size_t nj=0; nj<=degree-ni; ++nj) {
                    for (size_t nk=0; nk<=degree-ni-nj; ++nk) {
                        if (ni+nj+nk==degree) {
                            ORSA_DEBUG("global translated N[%i][%i][%i] = %g",ni,nj,nk,translated_N[ni][nj][nk]);
                        }
                    }
                }
            }
        }
    }
    {
        const double radiusCorrectionRatio = plateModelR0/gravityDataR0;
        norm_C.resize(SH_degree+1);
        norm_S.resize(SH_degree+1);
        for (size_t l=0; l<=SH_degree; ++l) {
            const double radiusCorrectionFactor = orsa::int_pow(radiusCorrectionRatio,l);
            norm_C[l].resize(l+1);
            norm_S[l].resize(l+1);
            for (size_t m=0; m<=l; ++m) {
                const orsa::triIndex_d C_tri_norm = orsa::conversionCoefficients_C_norm(l,m);
                const orsa::triIndex_d S_tri_norm = orsa::conversionCoefficients_S_norm(l,m);
                for (size_t ni=0; ni<=l; ++ni) {
                    for (size_t nj=0; nj<=l-ni; ++nj) {
                        for (size_t nk=0; nk<=l-ni-nj; ++nk) {
                            norm_C[l][m] +=
                                radiusCorrectionFactor *
                                C_tri_norm[ni][nj][nk] *
                                translated_N[ni][nj][nk] / translated_N[0][0][0];  
                            if (m != 0) {
                                norm_S[l][m] +=
                                    radiusCorrectionFactor *
                                    S_tri_norm[ni][nj][nk] *
                                    translated_N[ni][nj][nk] / translated_N[0][0][0];  
                            } else {
                                norm_S[l][m] = 0.0;
                            }
                        }
                    }
                }
            }
        }
        // debug output
        ORSA_DEBUG("reference R0: %g [km]",orsa::FromUnits(gravityDataR0,orsa::Unit::KM,-1));
        for (size_t l=0; l<=SH_degree; ++l) {
            for (size_t m=0; m<=l; ++m) {
                ORSA_DEBUG("norm_C[%i][%i] = %Fg",l,m,norm_C[l][m].get_mpf_t());
                if (m != 0) ORSA_DEBUG("norm_S[%i][%i] = %Fg",l,m,norm_S[l][m].get_mpf_t());
            }
        }
    }
}

#endif // CCMD2SH_H
