#ifndef CCMD2SH_H
#define CCMD2SH_H

#include "CubicChebyshevMassDistribution.h"
#include "simplex.h"
#include "SH2ijk.h"
#include <orsa/paulMoment.h>

// note on CM: it is used as input if set, or it is computed using CCMD if unset
template <typename T>
void CCMD2SH(orsa::Cache<orsa::Vector> & CM,
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

    const bool verbose=false;
    
    double mass_cT = 0.0;
    double i1d = 0.0;
    {
        orsa::Vector CM_sum_vector(0,0,0);
        double CM_sum_mass = 0.0;
        //
        double CMx_over_plateModelR0 = 0.0;
        double CMy_over_plateModelR0 = 0.0;
        double CMz_over_plateModelR0 = 0.0;
        {
            // double i1d = 0.0;
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
            if (i1d != 0.0) {
                mass_cT = i1d*orsa::cube(plateModelR0)*CCMD->densityScale;    
                ORSA_DEBUG("mass_cT: %g",mass_cT);
                CMx_over_plateModelR0 = iXd / i1d;
                CMy_over_plateModelR0 = iYd / i1d;
                CMz_over_plateModelR0 = iZd / i1d;
                //
                CM_sum_vector +=
                    mass_cT *
                    orsa::Vector(CMx_over_plateModelR0*plateModelR0,
                                 CMy_over_plateModelR0*plateModelR0,
                                 CMz_over_plateModelR0*plateModelR0);
                CM_sum_mass += mass_cT;
            }
        }
        if (CCMD->layerData.get() != 0) {
            {
                const LayerData::EllipsoidLayerVectorType & elv = CCMD->layerData->ellipsoidLayerVector;
                ORSA_DEBUG("elv.size(): %i",elv.size());
                for (size_t k=0; k<elv.size(); ++k) {
                    const double elv_excessMass = elv[k]->volume()*elv[k]->excessDensity;
                    CM_sum_vector += elv_excessMass * elv[k]->v0;
                    CM_sum_mass   += elv_excessMass;
                }
            }
            {
                const LayerData::SHLayerVectorType & shlv = CCMD->layerData->shLayerVector;
                ORSA_DEBUG("shlv.size(): %i",shlv.size());
                for (size_t k=0; k<shlv.size(); ++k) {
                    const double shlv_excessMass = shlv[k]->volume()*shlv[k]->excessDensity;
                    CM_sum_vector += shlv_excessMass * shlv[k]->v0;
                    CM_sum_mass   += shlv_excessMass;
                }
            }
        }
        // CM is also an output variable
        if (!CM.isSet()) CM = CM_sum_vector / CM_sum_mass;
        print(CM);
    }
    
    const double CMx_over_plateModelR0 = (*CM).getX()/plateModelR0;
    const double CMy_over_plateModelR0 = (*CM).getY()/plateModelR0;
    const double CMz_over_plateModelR0 = (*CM).getZ()/plateModelR0;
    
    // determine total mass, which is the sum of the contributions form Chebyshev (mass_cT) and Layers (below)
    double mass_layer = 0.0;
    if (CCMD->layerData.get() != 0) {
        {
            const LayerData::EllipsoidLayerVectorType & elv = CCMD->layerData->ellipsoidLayerVector;
            for (size_t k=0; k<elv.size(); ++k) {
                mass_layer += elv[k]->excessDensity*elv[k]->volume();
                
                ORSA_DEBUG("elv[%i]->excessDensity: %g   elv[%i]->volume: %g",
                           k,elv[k]->excessDensity,
                           k,elv[k]->volume());
            }
        }
        {
            const LayerData::SHLayerVectorType & shlv = CCMD->layerData->shLayerVector;
            for (size_t k=0; k<shlv.size(); ++k) {
                mass_layer += shlv[k]->excessDensity*shlv[k]->volume();
                
                ORSA_DEBUG("shlv[%i]->excessDensity: %g   shlv[%i]->volume: %g",
                           k,shlv[k]->excessDensity,
                           k,shlv[k]->volume());
            }
        }
    }
    const double totalMass = mass_cT + mass_layer;
    
    ORSA_DEBUG("totalMass: %g",totalMass);
    
    norm_C.resize(SH_degree+1);
    norm_S.resize(SH_degree+1);
    for (size_t l=0; l<=SH_degree; ++l) {
        norm_C[l].resize(l+1);
        norm_S[l].resize(l+1);
        for (size_t m=0; m<=l; ++m) {
            norm_C[l][m] = 0.0;
            norm_S[l][m] = 0.0;
        }
    }
    
    if (totalMass == 0.0) {
        ORSA_DEBUG("problem: totalMass = %g",totalMass);
        return;
    }
    
    const double cT_massFactor = mass_cT / totalMass;
    
    ORSA_DEBUG("cT_massFactor: %g",cT_massFactor);
    
    if (cT_massFactor != 0.0) for (size_t l=0; l<=SH_degree; ++l) {
        const double radiusCorrectionFactor = orsa::int_pow(radiusCorrectionRatio,l);
        for (size_t m=0; m<=l; ++m) {
            const orsa::triIndex_mpq C_tri_integral = orsa::conversionCoefficients_C_integral(l,m);
            const orsa::triIndex_d   C_tri_norm     = orsa::conversionCoefficients_C_norm(l,m);
            // norm_C[l][m] = 0.0;
            const orsa::triIndex_mpq S_tri_integral = orsa::conversionCoefficients_S_integral(l,m);
            const orsa::triIndex_d   S_tri_norm     = orsa::conversionCoefficients_S_norm(l,m);
            // norm_S[l][m] = 0.0;
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
            
            norm_C[l][m] *= radiusCorrectionFactor;
            norm_S[l][m] *= radiusCorrectionFactor;
            
            norm_C[l][m] /= i1d;
            norm_S[l][m] /= i1d;
            
            norm_C[l][m] *= cT_massFactor;
            norm_S[l][m] *= cT_massFactor;
            
            ORSA_DEBUG("norm_C[%i][%i] = %Fg  [cT only]",l,m,norm_C[l][m].get_mpf_t());
            if (m != 0) ORSA_DEBUG("norm_S[%i][%i] = %Fg  [cT only]",l,m,norm_S[l][m].get_mpf_t());
            
        }
    }
    
    if (CCMD->layerData.get() != 0) {
        
        // add layerData contribution
        
        // now, layer by layer
        {
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
                
                const double elv_excessMass = elv[k]->volume()*elv[k]->excessDensity;
                ORSA_DEBUG("elv_excessMass[%i]: %g",k,elv_excessMass);
                const double elv_massFactor = elv_excessMass / totalMass;
                ORSA_DEBUG("elv_massFactor[%i]: %g",k,elv_massFactor);
                
                // scale contribution and add to total Clm and Slm
                for (size_t l=0; l<=SH_degree; ++l) {
                    for (size_t m=0; m<=l; ++m) {
                        // scale by mass
                        elv_norm_C[l][m] *= elv_massFactor;
                        elv_norm_S[l][m] *= elv_massFactor;
                        
                        ORSA_DEBUG("elv[%i]_norm_C[%i][%i] = %Fg",k,l,m,elv_norm_C[l][m].get_mpf_t());
                        if (m != 0) ORSA_DEBUG("elv[%i]_norm_S[%i][%i] = %Fg",k,l,m,elv_norm_S[l][m].get_mpf_t());
                        
                        norm_C[l][m] += elv_norm_C[l][m];
                        norm_S[l][m] += elv_norm_S[l][m];
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
                                                                            SQLiteDBFileName);
                
                // const size_t layer_degree = shlv[k]->norm_A.size()-1;
                
                std::vector< std::vector<mpf_class> > shlv_norm_C;
                std::vector< std::vector<mpf_class> > shlv_norm_S;
                shlv_norm_C.resize(SH_degree+1);
                shlv_norm_S.resize(SH_degree+1);
                for (size_t l=0; l<=SH_degree; ++l) {
                    shlv_norm_C[l].resize(l+1);
                    shlv_norm_S[l].resize(l+1);
                     for (size_t m=0; m<=l; ++m) {
                         shlv_norm_C[l][m] = 0.0;
                         shlv_norm_S[l][m] = 0.0;
                     }
                }
                
                const double CMx_over_plateModelR0 = shlv[k]->v0.getX()/plateModelR0;
                const double CMy_over_plateModelR0 = shlv[k]->v0.getY()/plateModelR0;
                const double CMz_over_plateModelR0 = shlv[k]->v0.getZ()/plateModelR0;
                
                for (size_t l=0; l<=SH_degree; ++l) {
                    const double radiusCorrectionFactor = orsa::int_pow(radiusCorrectionRatio,l);
                    for (size_t m=0; m<=l; ++m) {
                        
                        // ORSA_DEBUG("l=%i   m=%i",l,m);
                        
                        const orsa::triIndex_mpq C_tri_integral = orsa::conversionCoefficients_C_integral(l,m);
                        const orsa::triIndex_d   C_tri_norm     = orsa::conversionCoefficients_C_norm(l,m);
                        
                        const orsa::triIndex_mpq S_tri_integral = orsa::conversionCoefficients_S_integral(l,m);
                        const orsa::triIndex_d   S_tri_norm     = orsa::conversionCoefficients_S_norm(l,m);
                        
                        /* const size_t z_C = (l==0) ?
                           mod_gravityData_index(gravityData.get(),"GM") :
                           mod_gravityData_index(gravityData.get(),orsaPDS::RadioScienceGravityData::keyC(l,m));
                           const size_t z_S = (m==0) ? 0 : mod_gravityData_index(gravityData.get(),orsaPDS::RadioScienceGravityData::keyS(l,m));
                        */
                        
                        // ni,nj,nk are the expansion of C_lm,S_lm in terms of N_ijk
                        for (size_t ni=0; ni<=l; ++ni) {
                            for (size_t nj=0; nj<=l-ni; ++nj) {
                                for (size_t nk=0; nk<=l-ni-nj; ++nk) {
                                    if ( (C_tri_integral[ni][nj][nk] == 0) && (S_tri_integral[ni][nj][nk] == 0) ) continue;
                                    // ti,tj,tk are the expansion of the density in terms of the cubic Chebyshev
                                    /* for (size_t running_T_degree=0; running_T_degree<=T_degree; ++running_T_degree) {
                                       for (size_t ti=0; ti<=T_degree; ++ti) {
                                       for (size_t tj=0; tj<=T_degree-ti; ++tj) {
                                       for (size_t tk=0; tk<=T_degree-ti-tj; ++tk) {
                                       if (ti+tj+tk != running_T_degree) continue;
                                       const std::vector<mpz_class> & cTi = orsa::ChebyshevTcoeff(ti);
                                       const std::vector<mpz_class> & cTj = orsa::ChebyshevTcoeff(tj);
                                       const std::vector<mpz_class> & cTk = orsa::ChebyshevTcoeff(tk);
                                       
                                       const size_t z_cT = CubicChebyshevMassDistribution::index(ti,tj,tk);
                                       
                                       double C2cT = 0.0;
                                       double S2cT = 0.0;
                                       
                                       // ci,cj,ck are the expansion of each Chebyshev polynomial in terms of powers of x,y,z
                                       for (size_t ci=0; ci<=ti; ++ci) {
                                       if (cTi[ci] == 0) continue;
                                       for (size_t cj=0; cj<=tj; ++cj) {
                                       if (cTj[cj] == 0) continue;
                                       for (size_t ck=0; ck<=tk; ++ck) {
                                       if (cTk[ck] == 0) continue;
                                    */
                                    // bi,bj,bk are the binomial expansion about the center of mass
                                    // this also introduces a power_sign
                                    for (size_t bi=0; bi<=ni; ++bi) {
                                        for (size_t bj=0; bj<=nj; ++bj) {
                                            for (size_t bk=0; bk<=nk; ++bk) {

                                                shlv_norm_C[l][m] +=
                                                    orsa::power_sign(bi+bj+bk) *
                                                    radiusCorrectionFactor *
                                                    C_tri_norm[ni][nj][nk] *
                                                    mpz_class(orsa::binomial(ni,bi) *
                                                              orsa::binomial(nj,bj) *
                                                              orsa::binomial(nk,bk)).get_d() *
                                                    orsa::int_pow(CMx_over_plateModelR0,bi) *
                                                    orsa::int_pow(CMy_over_plateModelR0,bj) *
                                                    orsa::int_pow(CMz_over_plateModelR0,bk) *
                                                    shi->getIntegral(ni-bi,nj-bj,nk-bk,verbose);
                                                
                                                shlv_norm_S[l][m] += 
                                                    orsa::power_sign(bi+bj+bk) *
                                                    radiusCorrectionFactor *
                                                    S_tri_norm[ni][nj][nk] *
                                                    mpz_class(orsa::binomial(ni,bi) *
                                                              orsa::binomial(nj,bj) *
                                                              orsa::binomial(nk,bk)).get_d() *
                                                    orsa::int_pow(CMx_over_plateModelR0,bi) *
                                                    orsa::int_pow(CMy_over_plateModelR0,bj) *
                                                    orsa::int_pow(CMz_over_plateModelR0,bk) *
                                                    shi->getIntegral(ni-bi,nj-bj,nk-bk,verbose);
                                                
                                                /* 
                                                   if (C_tri_integral[ni][nj][nk] != 0) {
                                                   C2cT +=
                                                   orsa::power_sign(bi+bj+bk) *
                                                   radiusCorrectionFactor *
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
                                                   S2cT +=
                                                   orsa::power_sign(bi+bj+bk) *
                                                   radiusCorrectionFactor *
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
                                                */
                                                
                                            }
                                        }
                                    }
                                    /* }
                                       }                                        
                                       }
                                       
                                       C2cT /= si->getIntegral(0,0,0);
                                       if (m!=0) S2cT /= si->getIntegral(0,0,0);
                                       
                                       if (l==0) {
                                       C2cT *= GM;
                                       }
                                       
                                       // saving, NOTE how we are adding terms here
                                       gsl_matrix_set(cT2sh,z_C,z_cT,gsl_matrix_get(cT2sh,z_C,z_cT)+C2cT);
                                       if (m!=0) gsl_matrix_set(cT2sh,z_S,z_cT,gsl_matrix_get(cT2sh,z_S,z_cT)+S2cT);
                                       }
                                       }
                                       }
                                       }
                                    */
                                }
                            }
                        }
                    }
                }
                
                const double shlv_excessMass = shlv[k]->volume()*shlv[k]->excessDensity;
                ORSA_DEBUG("shlv_excessMass[%i]: %g",k,shlv_excessMass);
                const double shlv_massFactor = shlv_excessMass / totalMass;
                ORSA_DEBUG("shlv_massFactor[%i]: %g",k,shlv_massFactor);
                
                for (size_t l=0; l<=SH_degree; ++l) {
                    for (size_t m=0; m<=l; ++m) {
                        shlv_norm_C[l][m] *= shlv_massFactor;
                        shlv_norm_S[l][m] *= shlv_massFactor;
                        
                        shlv_norm_C[l][m] /= shi->getIntegral(0,0,0);
                        shlv_norm_S[l][m] /= shi->getIntegral(0,0,0);
                        
                        norm_C[l][m] += shlv_norm_C[l][m];
                        norm_S[l][m] += shlv_norm_S[l][m];
                        
                        ORSA_DEBUG("shlv_norm_C[%i][%i] = %Fg",l,m,shlv_norm_C[l][m].get_mpf_t());
                        if (m != 0) ORSA_DEBUG("shlv_norm_S[%i][%i] = %Fg",l,m,shlv_norm_S[l][m].get_mpf_t());
                    }
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
