#ifndef CCMD2IJK_H
#define CCMD2IJK_H

#include "CubicChebyshevMassDistribution.h"
#include "simplex.h"
#include "SH2ijk.h"
#include "global_SH_epsrel.h"

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
                                        
                                        ORSA_DEBUG("N[%i][%i][%i] += %12.6g x %12.6g     t: %i %i %i   c: %i %i %i",
                                                   ni,nj,nk,
                                                   baseFactor,
                                                   si->getIntegral(ci+ni,cj+nj,ck+nk),
                                                   ti,tj,tk,
                                                   ci,cj,ck);
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
                for (size_t ni=0; ni<=degree; ++ni) {
                    for (size_t nj=0; nj<=degree-ni; ++nj) {
                        for (size_t nk=0; nk<=degree-ni-nj; ++nk) {
                            // bi,bj,bk are the binomial expansion of the v0 translation
                            // this also introduces a power_sign
                            for (size_t bi=0; bi<=ni; ++bi) {
                                for (size_t bj=0; bj<=nj; ++bj) {
                                    for (size_t bk=0; bk<=nk; ++bk) {
                                        
                                        /* if ((ni-bi)%2==1) continue;
                                           if ((nj-bj)%2==1) continue;
                                           if ((nk-bk)%2==1) continue;
                                        */
                                        
                                        const mpz_class factor_i   = EllipsoidExpansion_product_utility((ni-bi)/2);
                                        const mpz_class factor_j   = EllipsoidExpansion_product_utility((nj-bj)/2);
                                        const mpz_class factor_k   = EllipsoidExpansion_product_utility((nk-bk)/2);
                                        const mpz_class factor_ijk = EllipsoidExpansion_product_utility((ni+nj+nk-bi-bj-bk)/2+2);
                                        //
                                        mpq_class factor(3*(factor_i*factor_j*factor_k),factor_ijk);
                                        factor.canonicalize();
                                        
                                        N[ni][nj][nk] +=
                                            elv[k]->excessDensity *
                                            orsa::power_sign(bi+bj+bk) *
                                            mpz_class(orsa::binomial(ni,bi) *
                                                      orsa::binomial(nj,bj) *
                                                      orsa::binomial(nk,bk)).get_d() *
                                            orsa::int_pow(elv[k]->v0.getX(),bi) *
                                            orsa::int_pow(elv[k]->v0.getY(),bj) *
                                            orsa::int_pow(elv[k]->v0.getZ(),bk) *
                                            factor.get_d() *
                                            orsa::int_pow(elv[k]->a/plateModelR0,ni-bi) *
                                            orsa::int_pow(elv[k]->b/plateModelR0,nj-bj) *
                                            orsa::int_pow(elv[k]->c/plateModelR0,nk-bk);
                                        
                                    }
                                }
                            }
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
                
                for (size_t ni=0; ni<=degree; ++ni) {
                    for (size_t nj=0; nj<=degree-ni; ++nj) {
                        for (size_t nk=0; nk<=degree-ni-nj; ++nk) {
                            // bi,bj,bk are the binomial expansion of the v0 translation
                            // this also introduces a power_sign
                            for (size_t bi=0; bi<=ni; ++bi) {
                                for (size_t bj=0; bj<=nj; ++bj) {
                                    for (size_t bk=0; bk<=nk; ++bk) {
                                        
                                        N[ni][nj][nk] +=
                                           shlv[k]->excessDensity *
                                            orsa::power_sign(bi+bj+bk) *
                                            mpz_class(orsa::binomial(ni,bi) *
                                                      orsa::binomial(nj,bj) *
                                                      orsa::binomial(nk,bk)).get_d() *
                                            orsa::int_pow(shlv[k]->v0.getX(),bi) *
                                            orsa::int_pow(shlv[k]->v0.getY(),bj) *
                                            orsa::int_pow(shlv[k]->v0.getZ(),bk) *
                                            shi->getIntegral(ni-bi,nj-bj,nk-bk,verbose);
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

#endif // CCMD2IJK_H
