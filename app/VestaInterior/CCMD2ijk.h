#ifndef CCMD2IJK_H
#define CCMD2IJK_H

#include "CubicChebyshevMassDistribution.h"
#include "simplex.h"
#include "SH2ijk.h"
#include "global_SH_epsrel.h"
#include "translate_ijk.h"

static mpz_class EllipsoidExpansion_product_utility(const unsigned int n) {
    mpz_class product = 1;
    for (unsigned int k=1; k<=n; ++k) {
        product *= (2*k-1);
    }
    return product;
}

#warning CHECK the signs of the binomial expansions!

template <typename T>
void CCMD2ijk(std::vector< std::vector< std::vector<double> > > & N,
              const size_t                                      & degree,
              const SimplexIntegration<T>                       * si,
              const CubicChebyshevMassDistribution              * CCMD,
              const double                                      & plateModelR0) {
    
    N.resize(degree+1);
    for (size_t ni=0; ni<=degree; ++ni) {
        N[ni].resize(degree+1-ni);
        for (size_t nj=0; nj<=degree-ni; ++nj) {
            N[ni][nj].resize(degree+1-ni-nj);
            for (size_t nk=0; nk<=degree-ni-nj; ++nk) {
                N[ni][nj][nk] = 0.0;
            }
        }
    }
    
    const size_t T_degree = CCMD->degree();
    const double volume = si->getIntegral(0,0,0)*orsa::cube(plateModelR0);
    
    ORSA_DEBUG("volume: %g",volume);
    
    const bool verbose=true;
    
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
                                CCMD->densityScale *
                                CCMD->coeff[ti][tj][tk] *
                                mpz_class(cTi[ci] * cTj[cj] * cTk[ck]).get_d();
                            for (size_t ni=0; ni<=degree; ++ni) {
                                for (size_t nj=0; nj<=degree-ni; ++nj) {
                                    for (size_t nk=0; nk<=degree-ni-nj; ++nk) {
                                        N[ni][nj][nk] += baseFactor * si->getIntegral(ci+ni,cj+nj,ck+nk);
                                        /* ORSA_DEBUG("N[%i][%i][%i] += %12.6g x %12.6g     t: %i %i %i   c: %i %i %i",
                                           ni,nj,nk,
                                           baseFactor,
                                           si->getIntegral(ci+ni,cj+nj,ck+nk),
                                           ti,tj,tk,
                                           ci,cj,ck);
                                        */
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    if (CCMD->layerData.get() != 0) {
        
        // add layerData contribution
        
        // now, layer by layer
        {
            const LayerData::EllipsoidLayerVectorType & elv = CCMD->layerData->ellipsoidLayerVector;
            for (size_t k=0; k<elv.size(); ++k) {

                std::vector< std::vector< std::vector<double> > > N_ell;
                N_ell.resize(degree+1);
                for (size_t ni=0; ni<=degree; ++ni) {
                    N_ell[ni].resize(degree+1-ni);
                    for (size_t nj=0; nj<=degree-ni; ++nj) {
                        N_ell[ni][nj].resize(degree+1-ni-nj);
                        for (size_t nk=0; nk<=degree-ni-nj; ++nk) {
                            if (ni%2==1) continue;
                            if (nj%2==1) continue;
                            if (nk%2==1) continue;
                            const mpz_class factor_i   = EllipsoidExpansion_product_utility(ni/2);
                            const mpz_class factor_j   = EllipsoidExpansion_product_utility(nj/2);
                            const mpz_class factor_k   = EllipsoidExpansion_product_utility(nk/2);
                            const mpz_class factor_ijk = EllipsoidExpansion_product_utility((ni+nj+nk)/2+2);
                            //
                            mpq_class triple_factor(factor_i*factor_j*factor_k,factor_ijk);
                            triple_factor.canonicalize();
                            //
                            N_ell[ni][nj][nk] +=
                                4*orsa::pi() *
                                elv[k]->excessDensity *
                                triple_factor.get_d() *
                                orsa::int_pow(elv[k]->a/plateModelR0,ni+1) *
                                orsa::int_pow(elv[k]->b/plateModelR0,nj+1) *
                                orsa::int_pow(elv[k]->c/plateModelR0,nk+1);
                        }
                    }
                }
                std::vector< std::vector< std::vector<double> > > translated_N_ell;
                translate(translated_N_ell,N_ell,elv[k]->v0);
                for (size_t ni=0; ni<=degree; ++ni) {
                    for (size_t nj=0; nj<=degree-ni; ++nj) {
                        for (size_t nk=0; nk<=degree-ni-nj; ++nk) {
                            N[ni][nj][nk] +=
                                translated_N_ell[ni][nj][nk];
                        }
                    }
                }
            }
        }
        
        {
            const LayerData::SHLayerVectorType & shlv = CCMD->layerData->shLayerVector;
            for (size_t k=0; k<shlv.size(); ++k) {
                const std::string SQLiteDBFileName = getSqliteDBFileName_SH(shlv[k]->MD5(),plateModelR0);
                osg::ref_ptr< SHIntegration<T> > shi = new SHIntegration<T>(shlv[k]->norm_A,
                                                                            shlv[k]->norm_B,
                                                                            plateModelR0,
                                                                            global_SH_epsrel,
                                                                            SQLiteDBFileName);
                std::vector< std::vector< std::vector<double> > > N_SH;
                N_SH.resize(degree+1);
                for (size_t ni=0; ni<=degree; ++ni) {
                    N_SH[ni].resize(degree+1-ni);
                    for (size_t nj=0; nj<=degree-ni; ++nj) {
                        N_SH[ni][nj].resize(degree+1-ni-nj);
                        for (size_t nk=0; nk<=degree-ni-nj; ++nk) {
                            N_SH[ni][nj][nk] =
                                shlv[k]->excessDensity *
                                shi->getIntegral(ni,nj,nk,verbose);
                        }
                    }
                }
                std::vector< std::vector< std::vector<double> > > translated_N_SH;
                translate(translated_N_SH,N_SH,shlv[k]->v0);
                for (size_t ni=0; ni<=degree; ++ni) {
                    for (size_t nj=0; nj<=degree-ni; ++nj) {
                        for (size_t nk=0; nk<=degree-ni-nj; ++nk) {
                            N[ni][nj][nk] +=
                                translated_N_SH[ni][nj][nk];
                        }
                    }
                }
            }
        }
    }
}

#endif // CCMD2IJK_H
