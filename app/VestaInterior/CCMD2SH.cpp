#include "CCMD2SH.h"

#include <orsa/paulMoment.h>

/*** CHOOSE ONE ***/
// typedef double simplex_T;
// typedef mpf_class simplex_T;
typedef dd_real simplex_T;
// typedef qd_real simplex_T;

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;

template <typename T>
void CCMD2SH(orsa::Vector & CM,
             std::vector< std::vector<mpf_class> > & norm_C,
             std::vector< std::vector<mpf_class> > & norm_S,
             const size_t                          & SH_degree,
             const SimplexIntegration<T>           * si,
             const CubicChebyshevMassDistribution::coeff & coeff,
             const double                          & plateModelR0,
             const double                          & gravityDataR0) {
    
    const size_t T_degree = CubicChebyshevMassDistribution::degree(coeff);
    const double radiusCorrectionRatio = plateModelR0/gravityDataR0;
    
    double i1d =0.0;
    double iXd =0.0;
    double iYd =0.0;
    double iZd =0.0;
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
                            const double baseFactor =
                                coeff[ti][tj][tk] *
                                mpz_class(cTi[ci] * cTj[cj] * cTk[ck]).get_d();
                            i1d  += baseFactor * si->getIntegral(ci,cj,ck);
                            iXd  += baseFactor * si->getIntegral(ci+1,cj,ck);
                            iYd  += baseFactor * si->getIntegral(ci,cj+1,ck);
                            iZd  += baseFactor * si->getIntegral(ci,cj,ck+1);
                        }
                    }
                }
            }
        }
    }
    const double CMx_over_plateModelR0 = iXd / i1d;
    const double CMy_over_plateModelR0 = iYd / i1d;
    const double CMz_over_plateModelR0 = iZd / i1d;
    
    // output
    CM = orsa::Vector(CMx_over_plateModelR0*plateModelR0,
                      CMy_over_plateModelR0*plateModelR0,
                      CMz_over_plateModelR0*plateModelR0);
    
    norm_C.resize(SH_degree+1);
    norm_S.resize(SH_degree+1);
    for (size_t l=0; l<=SH_degree; ++l) {
        norm_C[l].resize(l+1);
        norm_S[l].resize(l+1);
    }
    
    for (size_t l=0; l<=SH_degree; ++l) {
        const double radiusCorrectionFactor = orsa::int_pow(radiusCorrectionRatio,l);
        for (size_t m=0; m<=l; ++m) {
            const orsa::triIndex_mpq C_tri_integral = orsa::conversionCoefficients_C_integral(l,m);
            const orsa::triIndex_d   C_tri_norm     = orsa::conversionCoefficients_C_norm(l,m);
            norm_C[l][m] = 0.0;
            const orsa::triIndex_mpq S_tri_integral = orsa::conversionCoefficients_S_integral(l,m);
            const orsa::triIndex_d   S_tri_norm     = orsa::conversionCoefficients_S_norm(l,m);
            norm_S[l][m] = 0.0;
            // ni,nj,nk are the expansion of C_lm,S_lm in terms of N_ijk
            for (size_t ni=0; ni<=l; ++ni) {
                for (size_t nj=0; nj<=l-ni; ++nj) {
                    for (size_t nk=0; nk<=l-ni-nj; ++nk) {
                        if ( (C_tri_integral[ni][nj][nk] == 0) && (S_tri_integral[ni][nj][nk] == 0) ) continue;
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
                                                // bi,bj,bk are the binomial expansion about the center of mass
                                                // this also introduces a power_sign
                                                for (size_t bi=0; bi<=ni; ++bi) {
                                                    for (size_t bj=0; bj<=nj; ++bj) {
                                                        for (size_t bk=0; bk<=nk; ++bk) {
                                                            
                                                            if (C_tri_integral[ni][nj][nk] != 0) {
                                                                norm_C[l][m] +=
                                                                    orsa::power_sign(bi+bj+bk) *
                                                                    coeff[ti][tj][tk] *
                                                                    C_tri_norm[ni][nj][nk] *
                                                                    mpz_class(orsa::binomial(ni,bi) *
                                                                              orsa::binomial(nj,bj) *
                                                                              orsa::binomial(nk,bk) *
                                                                              cTi[ci] * cTj[cj] * cTk[ck]).get_d() *
                                                                    orsa::int_pow(CMx_over_plateModelR0,bi) *
                                                                    orsa::int_pow(CMy_over_plateModelR0,bj) *
                                                                    orsa::int_pow(CMz_over_plateModelR0,bk) *
                                                                    si->getIntegral(ni-bi+ci,nj-bj+cj,nk-bk+ck);
                                                            }
                                                            
                                                            if (S_tri_integral[ni][nj][nk] != 0) {
                                                                norm_S[l][m] +=
                                                                    orsa::power_sign(bi+bj+bk) *
                                                                    coeff[ti][tj][tk] *
                                                                    S_tri_norm[ni][nj][nk] *
                                                                    mpz_class(orsa::binomial(ni,bi) *
                                                                              orsa::binomial(nj,bj) *
                                                                              orsa::binomial(nk,bk) *
                                                                              cTi[ci] * cTj[cj] * cTk[ck]).get_d() *
                                                                    orsa::int_pow(CMx_over_plateModelR0,bi) *
                                                                    orsa::int_pow(CMy_over_plateModelR0,bj) *
                                                                    orsa::int_pow(CMz_over_plateModelR0,bk) *
                                                                    si->getIntegral(ni-bi+ci,nj-bj+cj,nk-bk+ck);
                                                            }                                                               
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            norm_C[l][m] /= si->getIntegral(0,0,0);
            norm_C[l][m] *= radiusCorrectionFactor;
            norm_S[l][m] /= si->getIntegral(0,0,0);
            norm_S[l][m] *= radiusCorrectionFactor;
            
            ORSA_DEBUG("norm_C[%i][%i] = %g",l,m,norm_C[l][m]);
            if (m != 0) ORSA_DEBUG("norm_S[%i][%i] = %g",l,m,norm_S[l][m]);
        }
    }
}
