#ifndef CCMD2SH_H
#define CCMD2SH_H

#include "CubicChebyshevMassDistribution.h"
#include "simplex.h"

#include <orsa/paulMoment.h>

template <typename T>
void CCMD2SH(orsa::Vector & CM, /* note: CM is just an output variable */
             std::vector< std::vector<mpf_class> > & norm_C,
             std::vector< std::vector<mpf_class> > & norm_S,
             const size_t                          & SH_degree,
             const SimplexIntegration<T>           * si,
             const CubicChebyshevMassDistribution  * CCMD,
             const double                          & plateModelR0,
             const double                          & gravityDataR0) {
    
    const size_t T_degree = CCMD->degree();
    const double radiusCorrectionRatio = plateModelR0/gravityDataR0;
    const double volume = si->getIntegral(0,0,0)*orsa::cube(plateModelR0);
    
    ORSA_DEBUG("volume: %g",volume);
    
    double i1d = 0.0;
    double iXd = 0.0;
    double iYd = 0.0;
    double iZd = 0.0;
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
                                CCMD->coeff[ti][tj][tk] *
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
    const double mass_cT = i1d*orsa::cube(plateModelR0)*CCMD->densityScale;    
    ORSA_DEBUG("mass_cT: %g",mass_cT);
    const double CMx_over_plateModelR0 = iXd / i1d;
    const double CMy_over_plateModelR0 = iYd / i1d;
    const double CMz_over_plateModelR0 = iZd / i1d;
    // init with this data, and then update with layerData
    orsa::Vector CM_sum_vector =
        mass_cT *
        orsa::Vector(CMx_over_plateModelR0*plateModelR0,
                     CMy_over_plateModelR0*plateModelR0,
                     CMz_over_plateModelR0*plateModelR0);
    double CM_sum_mass = mass_cT;
    if (CCMD->layerData.get() != 0) {
        const LayerData::EllipsoidLayerVectorType & elv = CCMD->layerData->ellipsoidLayerVector;
        for (size_t k=0; k<elv.size(); ++k) {
            const double elv_excessMass = elv[k]->volume*elv[k]->excessDensity;
            CM_sum_vector += elv_excessMass * elv[k]->v0;
            CM_sum_mass   += elv_excessMass;
        }
    }
    // CM is also an output variable
    CM = CM_sum_vector / CM_sum_mass;
    print(CM);
    
    // determine total mass, which is the sum of the contributions form Chebyshev (mass_cT) and Layers (below)
    double mass_layer = 0.0;
    if (CCMD->layerData.get() != 0) {
        // mass_layer += volume*CCMD->layerData->baseDensity;
        const LayerData::EllipsoidLayerVectorType & elv = CCMD->layerData->ellipsoidLayerVector;
        for (size_t k=0; k<elv.size(); ++k) {
            mass_layer += elv[k]->excessDensity*elv[k]->volume;
            
            ORSA_DEBUG("elv[%i]->excessDensity: %g   elv[%i]->volume: %g",
                       k,elv[k]->excessDensity,
                       k,elv[k]->volume);
        }
    }
    const double totalMass = mass_cT + mass_layer;
    
    ORSA_DEBUG("totalMass: %g",totalMass);
    
    norm_C.resize(SH_degree+1);
    norm_S.resize(SH_degree+1);
    for (size_t l=0; l<=SH_degree; ++l) {
        norm_C[l].resize(l+1);
        norm_S[l].resize(l+1);
    }

    const double cT_massFactor = mass_cT / totalMass;
    
    ORSA_DEBUG("cT_massFactor: %g",cT_massFactor);
    
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
                                                                    CCMD->coeff[ti][tj][tk] *
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
                                                                    CCMD->coeff[ti][tj][tk] *
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
            
            norm_C[l][m] *= cT_massFactor;
            norm_S[l][m] *= cT_massFactor;
            
        }
    }
    
    if (CCMD->layerData.get() != 0) {
        
        // add layerData contribution
        
        // first, the baseDensity (bD) contribution, which is proportional to the shape
        /* if (CCMD->layerData->baseDensity != 0.0) {
           #warning use epsilon instead of 0.0 ?
           orsa::Vector bD_CM;
           std::vector< std::vector<mpf_class> > bD_norm_C;
           std::vector< std::vector<mpf_class> > bD_norm_S;
           
           CubicChebyshevMassDistribution::CoefficientType bD_coeff;
           CubicChebyshevMassDistribution::resize(bD_coeff,0);
           bD_coeff[0][0][0] = 1.0;
           
           osg::ref_ptr<CubicChebyshevMassDistribution> bD_CCMD =
           new CubicChebyshevMassDistribution(bD_coeff,
           CCMD->layerData->baseDensity,
           plateModelR0,
           0);
           
           CCMD2SH(bD_CM,
           bD_norm_C,
           bD_norm_S,
           SH_degree,
           si,
           bD_CCMD.get(),
           plateModelR0,
           gravityDataR0);
           
           #warning CHECK THIS
           const double bD_excessMass = volume*CCMD->layerData->baseDensity;
           const double bD_massFactor = bD_excessMass / totalMass;
           ORSA_DEBUG("bD_massFactor: %g",bD_massFactor);
           
           // scale contribution and add to total Clm and Slm
           for (size_t l=0; l<=SH_degree; ++l) {
           for (size_t m=0; m<=l; ++m) {
           // scale by mass
           bD_norm_C[l][m] *= bD_massFactor;
           bD_norm_S[l][m] *= bD_massFactor;
           
           ORSA_DEBUG("bD_norm_C[%i][%i] = %Fg",l,m,bD_norm_C[l][m].get_mpf_t());
           if (m != 0) ORSA_DEBUG("bD_norm_S[%i][%i] = %Fg",l,m,bD_norm_S[l][m].get_mpf_t());
           
           norm_C[l][m] += bD_norm_C[l][m];
           norm_S[l][m] += bD_norm_S[l][m];
           }
           }
           }
        */
        
        // now, layer by layer
        osg::ref_ptr<orsa::PaulMoment> elv_pm = new orsa::PaulMoment(SH_degree);
        osg::ref_ptr<orsa::PaulMoment> elv_translated_pm = new orsa::PaulMoment(SH_degree);
        std::vector< std::vector<mpf_class> > elv_C;
        std::vector< std::vector<mpf_class> > elv_S;
        std::vector< std::vector<mpf_class> > elv_norm_C;
        std::vector< std::vector<mpf_class> > elv_norm_S;
        std::vector<mpf_class> elv_J;
        const LayerData::EllipsoidLayerVectorType & elv = CCMD->layerData->ellipsoidLayerVector;
        for (size_t k=0; k<elv.size(); ++k) {
            orsa::EllipsoidExpansion(elv_pm.get(),
                                     elv[k]->a,
                                     elv[k]->b,
                                     elv[k]->c);
            
            if (!translate(elv_translated_pm.get(),
                           elv_pm.get(),
                           elv[k]->v0)) {
                ORSA_DEBUG("problems...");
                exit(0);
            }
            
            orsa::convert(elv_C,
                          elv_S,
                          elv_norm_C,
                          elv_norm_S,
                          elv_J,
                          elv_translated_pm.get(),
                          gravityDataR0);
            
            const double elv_excessMass = elv[k]->volume*elv[k]->excessDensity;
            ORSA_DEBUG("elv_excessMass[%i]: %g",k,elv_excessMass);
            const double elv_massFactor = elv_excessMass / totalMass;
            ORSA_DEBUG("elv_massFactor[%i]: %g",k,elv_massFactor);
            
            // scale contribution and add to total Clm and Slm
            for (size_t l=0; l<=SH_degree; ++l) {
                for (size_t m=0; m<=l; ++m) {
                    // scale by mass
                    elv_norm_C[l][m] *= elv_massFactor;
                    elv_norm_S[l][m] *= elv_massFactor;
                    
                    /* ORSA_DEBUG("elv_norm_C[%i][%i] = %Fg",l,m,elv_norm_C[l][m].get_mpf_t());
                       if (m != 0) ORSA_DEBUG("elv_norm_S[%i][%i] = %Fg",l,m,elv_norm_S[l][m].get_mpf_t());
                    */
                    
                    norm_C[l][m] += elv_norm_C[l][m];
                    norm_S[l][m] += elv_norm_S[l][m];
                }
            }
        }
    }
    
    // debug output
    for (size_t l=0; l<=SH_degree; ++l) {
        for (size_t m=0; m<=l; ++m) {
            ORSA_DEBUG("norm_C[%i][%i] = %Fg",l,m,norm_C[l][m].get_mpf_t());
            if (m != 0) ORSA_DEBUG("norm_S[%i][%i] = %Fg",l,m,norm_S[l][m].get_mpf_t());
        }
    }
}

#endif // CCMD2SH_H
