#include "MovingLayers.h"

#include <orsa/chebyshev.h>
#include <orsa/statistic.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>

// #include "vesta.h"
#include "gaskell.h"

#include "simplex.h"

#include "CCMD2SH.h"

/*** CHOOSE ONE ***/
// typedef double simplex_T;
// typedef mpf_class simplex_T;
typedef dd_real simplex_T;
// typedef qd_real simplex_T;

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;

/*******/

// modified versions of RadioScienceGravityData calls, to include C10,C11,S11
unsigned int mod_gravityData_index(const orsaPDS::RadioScienceGravityData * gravityData,
                                   const QString & key) {
    unsigned int index;
    if (key == orsaPDS::RadioScienceGravityData::keyC(1,0)) {
        index = 1;
    } else if (key == orsaPDS::RadioScienceGravityData::keyC(1,1)) {
        index = 2;
    } else if (key == orsaPDS::RadioScienceGravityData::keyS(1,1)) {
        index = 3;
    } else {
        index = gravityData->index(key);
        if (index != 0) index += 3;
    }
    return index;
}

QString mod_gravityData_key(const orsaPDS::RadioScienceGravityData * gravityData,
                            const unsigned int & index) {
    if (index == 0) {
        return gravityData->key(index);
    } else if (index == 1) {
        return orsaPDS::RadioScienceGravityData::keyC(1,0);        
    } else if (index == 2) {
        return orsaPDS::RadioScienceGravityData::keyC(1,1);        
    } else if (index == 3) {
        return orsaPDS::RadioScienceGravityData::keyS(1,1);        
    } else {
        return gravityData->key(index-3);
    }
}
//
double mod_gravityData_getCoeff(const orsaPDS::RadioScienceGravityData * gravityData,
                                const QString & key) {
    double coeff;
    if ( (key == orsaPDS::RadioScienceGravityData::keyC(1,0)) ||  
         (key == orsaPDS::RadioScienceGravityData::keyC(1,1)) ||
         (key == orsaPDS::RadioScienceGravityData::keyS(1,1)) ) {
        coeff = 0.0;
    } else {
        coeff = gravityData->getCoeff(key);
    }
    return coeff;
}
//
/* double mod_gravityData_getCovar(const QString & key1, const QString & key2) {
   double covar;
   ORSA_DEBUG("complete here...");
   }
*/
//
unsigned int mod_gravityData_numberOfCoefficients(const orsaPDS::RadioScienceGravityData * gravityData) {
    return gravityData->numberOfCoefficients+3;
}

gsl_vector * mod_gravityData_getCoefficientVector(const orsaPDS::RadioScienceGravityData * gravityData) {
    gsl_vector * mu = gravityData->getCoefficientVector();
    gsl_vector * mod_mu = gsl_vector_alloc(mod_gravityData_numberOfCoefficients(gravityData));
    for (unsigned int k=0; k<mod_gravityData_numberOfCoefficients(gravityData); ++k) {
        if (k==0) {
            gsl_vector_set(mod_mu,k,gsl_vector_get(mu,k));
        } else if ( (k==1) || (k==2) || (k==3) ) {
            gsl_vector_set(mod_mu,k,0.0);
        } else {
            gsl_vector_set(mod_mu,k,gsl_vector_get(mu,k-3));
        }
    }
    return mod_mu;
}

gsl_matrix * mod_gravityData_getCovarianceMatrix(const orsaPDS::RadioScienceGravityData * gravityData) {
    gsl_matrix * covm = gravityData->getCovarianceMatrix();
    gsl_matrix * mod_covm = gsl_matrix_alloc(mod_gravityData_numberOfCoefficients(gravityData),
                                             mod_gravityData_numberOfCoefficients(gravityData));
    for (unsigned int l=0; l<mod_gravityData_numberOfCoefficients(gravityData); ++l) {
        for (unsigned int m=0; m<mod_gravityData_numberOfCoefficients(gravityData); ++m) { 
            if ((l==0) && (m==0)) {
                gsl_matrix_set(mod_covm,l,m,gsl_matrix_get(covm,l,m));
            } else if ((l==1) || (l==2) || (l==3) || (m==1) || (m==2) || (m==3)) {
                gsl_matrix_set(mod_covm,l,m,0.0);
            } else if ((l==0) && (m!=0)) {
                gsl_matrix_set(mod_covm,l,m,gsl_matrix_get(covm,l,m-3));
            } else if ((l!=0) && (m==0)) {
                gsl_matrix_set(mod_covm,l,m,gsl_matrix_get(covm,l-3,m));
            } else {
                gsl_matrix_set(mod_covm,l,m,gsl_matrix_get(covm,l-3,m-3));
            }
        }
    }
    return mod_covm;   
}

gsl_matrix * mod_gravityData_getInverseCovarianceMatrix(const orsaPDS::RadioScienceGravityData * gravityData) {
    gsl_matrix * inv_covm = gravityData->getInverseCovarianceMatrix();
    gsl_matrix * mod_inv_covm = gsl_matrix_alloc(mod_gravityData_numberOfCoefficients(gravityData),
                                             mod_gravityData_numberOfCoefficients(gravityData));
    for (unsigned int l=0; l<mod_gravityData_numberOfCoefficients(gravityData); ++l) {
        for (unsigned int m=0; m<mod_gravityData_numberOfCoefficients(gravityData); ++m) { 
            if ((l==0) && (m==0)) {
                gsl_matrix_set(mod_inv_covm,l,m,gsl_matrix_get(inv_covm,l,m));
            } else if ((l==1) || (l==2) || (l==3) || (m==1) || (m==2) || (m==3)) {
                gsl_matrix_set(mod_inv_covm,l,m,0.0);
            } else if ((l==0) && (m!=0)) {
                gsl_matrix_set(mod_inv_covm,l,m,gsl_matrix_get(inv_covm,l,m-3));
            } else if ((l!=0) && (m==0)) {
                gsl_matrix_set(mod_inv_covm,l,m,gsl_matrix_get(inv_covm,l-3,m));
            } else {
                gsl_matrix_set(mod_inv_covm,l,m,gsl_matrix_get(inv_covm,l-3,m-3));
            }
        }
    }
    return mod_inv_covm;   
}

/**********/

int main(int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    ORSA_DEBUG("PID: %i",getpid());
    
    // QD
    unsigned int oldcw;
    fpu_fix_start(&oldcw);
    
    //ORSA_DEBUG("current mpf precision: %i",mpf_get_default_prec());
    mpf_set_default_prec(128);
    // mpf_set_default_prec(256);
    // mpf_set_default_prec(512);
    // ORSA_DEBUG("updated mpf precision: %i",mpf_get_default_prec());
    
    if (argc != 14) {
        printf("Usage: %s <RadioScienceGravityFile> <plate-model-file> <plate-model-R0_km> <gravity-degree> <polynomial-degree> <CM-x_km> <CM-y_km> <CM-z_km> <CM-sigma-x_km> <CM-sigma-y_km> <CM-sigma-z_km> <num-sample-points> <CCMDF-input-file>\n",argv[0]);
        exit(0);
    }   
    
    const std::string radioScienceGravityFile = argv[1];
    const std::string plateModelFile = argv[2];
    const double plateModelR0 = orsa::FromUnits(atof(argv[3]),orsa::Unit::KM);
    const int gravityDegree = atoi(argv[4]);
    const int polynomialDegree = atoi(argv[5]);
    const double CM_x = orsa::FromUnits(atof(argv[6]),orsa::Unit::KM);
    const double CM_y = orsa::FromUnits(atof(argv[7]),orsa::Unit::KM);
    const double CM_z = orsa::FromUnits(atof(argv[8]),orsa::Unit::KM);
    const double CM_sx = orsa::FromUnits(atof(argv[9]),orsa::Unit::KM);
    const double CM_sy = orsa::FromUnits(atof(argv[10]),orsa::Unit::KM);
    const double CM_sz = orsa::FromUnits(atof(argv[11]),orsa::Unit::KM);
    const int numSamplePoints = atoi(argv[12]);
    const bool have_CCMDF_file = (argc == 14);
    const std::string CCMDF_filename = (argc == 14) ? argv[13] : "";
    
    // safer over NFS
    sqlite3_vfs_register(sqlite3_vfs_find("unix-dotfile"), 1);
    
    const std::string SQLiteDBFileName = getSqliteDBFileName_simplex(plateModelFile,plateModelR0);
    
    osg::ref_ptr<orsaPDS::RadioScienceGravityData> gravityData = new orsaPDS::RadioScienceGravityData;
    orsaPDS::RadioScienceGravityFile::read(gravityData.get(),radioScienceGravityFile,512,1518);
    
    osg::ref_ptr<GaskellPlateModel> shapeModel = new GaskellPlateModel;
    if (!shapeModel->read(plateModelFile)) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    osg::ref_ptr<SimplexIntegration<simplex_T> > si = new SimplexIntegration<simplex_T>(shapeModel.get(), plateModelR0, SQLiteDBFileName);
    si->reserve(polynomialDegree);
    
    const size_t SH_degree = gravityDegree; // shperical harmonics degree
    const size_t  T_degree = polynomialDegree; // chebyshev polynomials degree
    
    const size_t  SH_size = (SH_degree+1)*(SH_degree+1);
    const size_t   T_size = CubicChebyshevMassDistribution::totalSize(T_degree);
    
    const double GM = gravityData->GM; 
    
    ORSA_DEBUG("SH_size: %d   T_size: %d   mod_gravityData_numberOfCoefficients: %d",SH_size,T_size,mod_gravityData_numberOfCoefficients(gravityData.get()));
    
    if (T_size <= SH_size) {
        ORSA_DEBUG("this method works only when the problem is under-determined, exiting");
        exit(0);
    }
    
    if (gravityDegree > gravityData->degree) {
        ORSA_DEBUG("requested gravity degree [%i] is larger than the RadioScienceGravityFile [%s] degree [%i]; exiting",
                   gravityDegree,radioScienceGravityFile.c_str(),gravityData->degree);
        exit(0);
    }
    
    osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution;
    if (have_CCMDF_file) {
        
        CubicChebyshevMassDistributionFile::DataContainer CCMDF;
        CubicChebyshevMassDistributionFile::read(CCMDF,CCMDF_filename);
        if (CCMDF.size() == 0) {
            ORSA_DEBUG("empty CCMDF file: [%s]",CCMDF_filename.c_str());
            exit(0);
        }
        if (CCMDF.size() > 1) {
            ORSA_DEBUG("CCMDF [%s] should contain only one set of coefficients.",CCMDF_filename.c_str());
        }
        
        massDistribution = CCMD(CCMDF[CCMDF.size()-1]);
    }

    {
        bool good_CCMD=true;
        if (!have_CCMDF_file) good_CCMD=false;
        if (massDistribution.get() == 0) good_CCMD=false;
        if ( (massDistribution.get() != 0) &&
             (massDistribution->layerData.get() == 0) ) good_CCMD=false;
        if ( (massDistribution.get() != 0) &&
             (massDistribution->layerData.get() != 0) &&
             (massDistribution->layerData->shLayerVector.size() == 0) ) good_CCMD=false;
        if (!good_CCMD) {
            ORSA_DEBUG("CCMDF with SHLayers required as input file");
            exit(0);
        }
    }
    
#warning check if there is any ROTATION between reference systems
    
#warning check code for HIGH degree, might have to rewrite linear algebra...
    
    const orsa::Vector sampled_CM(CM_x+orsa::GlobalRNG::instance()->rng()->gsl_ran_gaussian(CM_sx),  
                                  CM_y+orsa::GlobalRNG::instance()->rng()->gsl_ran_gaussian(CM_sy),
                                  CM_z+orsa::GlobalRNG::instance()->rng()->gsl_ran_gaussian(CM_sz));
    
    const double CMx_over_plateModelR0 = sampled_CM.getX()/plateModelR0;
    const double CMy_over_plateModelR0 = sampled_CM.getY()/plateModelR0;
    const double CMz_over_plateModelR0 = sampled_CM.getZ()/plateModelR0;
    
    const double volume = si->getIntegral(0,0,0)*orsa::cube(plateModelR0);
    
    const double bulkDensity = GM/orsa::Unit::G()/volume;
    
    ORSA_DEBUG("bulkDensity: %g [g/cm^3]",orsa::FromUnits(orsa::FromUnits(bulkDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
    
    const double radiusCorrectionRatio = plateModelR0/gravityData->R0;
    
    gsl_matrix * cT2sh = gsl_matrix_calloc(SH_size,T_size);
    
    for (size_t l=0; l<=(size_t)gravityDegree; ++l) {
        // #warning should try to print degree 1 terms, just as a check 
        // #warning SKIP l=1 in production!!!!!!!!!!!!!!!!
        // if (l==1) continue;
        const double radiusCorrectionFactor = orsa::int_pow(radiusCorrectionRatio,l);
        for (size_t m=0; m<=l; ++m) {
            
            ORSA_DEBUG("l=%i   m=%i",l,m);
            
            const orsa::triIndex_mpq C_tri_integral = orsa::conversionCoefficients_C_integral(l,m);
            const orsa::triIndex_d   C_tri_norm     = orsa::conversionCoefficients_C_norm(l,m);
            
            const orsa::triIndex_mpq S_tri_integral = orsa::conversionCoefficients_S_integral(l,m);
            const orsa::triIndex_d   S_tri_norm     = orsa::conversionCoefficients_S_norm(l,m);
            
            const size_t z_C = (l==0) ?
                mod_gravityData_index(gravityData.get(),"GM") :
                mod_gravityData_index(gravityData.get(),orsaPDS::RadioScienceGravityData::keyC(l,m));
            const size_t z_S = (m==0) ? 0 : mod_gravityData_index(gravityData.get(),orsaPDS::RadioScienceGravityData::keyS(l,m));
            
            // ni,nj,nk are the expansion of C_lm,S_lm in terms of N_ijk
            for (size_t ni=0; ni<=l; ++ni) {
                for (size_t nj=0; nj<=l-ni; ++nj) {
                    for (size_t nk=0; nk<=l-ni-nj; ++nk) {
                        if ( (C_tri_integral[ni][nj][nk] == 0) && (S_tri_integral[ni][nj][nk] == 0) ) continue;
                        // ti,tj,tk are the expansion of the density in terms of the cubic Chebyshev
                        for (size_t running_T_degree=0; running_T_degree<=T_degree; ++running_T_degree) {
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
                                                    // bi,bj,bk are the binomial expansion about the center of mass
                                                    // this also introduces a power_sign
                                                    for (size_t bi=0; bi<=ni; ++bi) {
                                                        for (size_t bj=0; bj<=nj; ++bj) {
                                                            for (size_t bk=0; bk<=nk; ++bk) {
                                                                
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
                                                            }
                                                        }
                                                    }
                                                }
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
                    }
                }
            }
        }
    }
    
    gsl_vector * pds_coeff    = mod_gravityData_getCoefficientVector(gravityData.get());
    gsl_matrix * pds_covm     = mod_gravityData_getCovarianceMatrix(gravityData.get());
    gsl_matrix * pds_inv_covm = mod_gravityData_getInverseCovarianceMatrix(gravityData.get());
    
    gsl_matrix_view     pds_covm_view = gsl_matrix_submatrix(    pds_covm, 0, 0, SH_size, SH_size);
    gsl_matrix_view pds_inv_covm_view = gsl_matrix_submatrix(pds_inv_covm, 0, 0, SH_size, SH_size);
    
    {
        
        // A x = b
        // A = cT2sh
        // x = A^T (A A^T)^(-1) b
        
        const size_t M = SH_size;
        const size_t N =  T_size;
        
        // A
        gsl_matrix * A = gsl_matrix_alloc(M,N);
        
        for (size_t j=0; j<M; ++j) {
            for (size_t k=0; k<N; ++k) {
                gsl_matrix_set(A,j,k,gsl_matrix_get(cT2sh,j,k));
            }
        }
        
        // A^T
        gsl_matrix * AT = gsl_matrix_alloc(N,M);
        
        for (size_t j=0; j<N; ++j) {
            for (size_t k=0; k<M; ++k) {
                gsl_matrix_set(AT,j,k,gsl_matrix_get(cT2sh,k,j));
            }
        }

        // QR decomposition of A^T, to find basis of null space
        
        gsl_matrix * QR = gsl_matrix_alloc(N,M);
        gsl_vector * tau = gsl_vector_alloc(std::min(M,N));
        
        gsl_matrix_memcpy(QR,AT);
        
        gsl_linalg_QR_decomp(QR,tau);

        gsl_matrix * Q = gsl_matrix_alloc(N,N);
        gsl_matrix * R = gsl_matrix_alloc(N,M);
        
        gsl_linalg_QR_unpack(QR,tau,Q,R);

        // null space basis
        gsl_vector * uK[N-M];
        for (size_t b=0; b<(N-M); ++b) {
            uK[b] = gsl_vector_alloc(N);
            for (size_t s=0; s<N; ++s) {
                gsl_vector_set(uK[b],s,gsl_matrix_get(Q,s,M+b));
                //
                if (0) {
                    size_t Tx,Ty,Tz;
                    CubicChebyshevMassDistribution::triIndex(Tx,Ty,Tz,s);
                    ORSA_DEBUG("uK[%03i][%03i] =%+12.6f (null space base vector for cT[%i][%i][%i])",b,s,gsl_vector_get(uK[b],s),Tx,Ty,Tz);
                }
            }
        }
        
        // (A A^T)
        gsl_matrix * A_AT = gsl_matrix_alloc(M,M);
        
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,AT,0.0,A_AT);
        
        /* for (size_t j=0; j<M; ++j) {
           for (size_t k=0; k<M; ++k) {
           ORSA_DEBUG("(A A^T)[%03i][%03i] = %+12.6g",j,k,gsl_matrix_get(A_AT,j,k));
           }
           }
        */
        
        // compute (A_AT)^(-1)
        gsl_matrix * inv_A_AT = gsl_matrix_alloc(M,M);
        {
            // first, find eigenvectors and eigenvalues of covariance matrix
            gsl_eigen_symmv_workspace * workspace = gsl_eigen_symmv_alloc(M);
            gsl_matrix * evec = gsl_matrix_alloc(M,M);
            gsl_vector * eval = gsl_vector_alloc(M);
            gsl_eigen_symmv(A_AT,eval,evec,workspace);
            // check if definite positive
            for (size_t k=0;k<M;++k) {
                if (gsl_vector_get(eval,k) <= 0.0) {
                    ORSA_DEBUG("problems... eval[%03i] = %g",k,gsl_vector_get(eval,k));
                    return 0;
                }
            }
            // inverse of eval matrix (diagonal with values 1/eigenval)
            gsl_matrix * inverse_eval_matrix = gsl_matrix_alloc(M,M);
            for (size_t i=0;i<M;++i) {
                for (size_t j=0;j<M;++j) {
                    if (i==j) {
                        gsl_matrix_set(inverse_eval_matrix,i,j,1.0/gsl_vector_get(eval,i));
                    } else {
                        gsl_matrix_set(inverse_eval_matrix,i,j,0.0);
                    }
                }
            }
            gsl_matrix * tmp_matrix = gsl_matrix_alloc(M,M);
            gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,inverse_eval_matrix,evec,0.0,tmp_matrix);
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,evec,tmp_matrix,0.0,inv_A_AT);
            
            // free all gsl allocated vars, except inv_A_AT
            gsl_eigen_symmv_free(workspace);    
            gsl_matrix_free(evec);
            gsl_vector_free(eval);
            // gsl_matrix_free(eval_matrix);
            gsl_matrix_free(inverse_eval_matrix);
            // skip inv_A_AT
            gsl_matrix_free(tmp_matrix);
        }

        // pseudo inverse of A = A^T (A A^T)^(-1)
        gsl_matrix * pseudoInvA = gsl_matrix_alloc(N,M);
        
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,AT,inv_A_AT,0.0,pseudoInvA);
        
        
        gsl_vector * sh = gsl_vector_calloc(M); // b  of A x = b
        gsl_vector * cT = gsl_vector_calloc(N); // x  of A x = b
        
        // sample from SH covariance matrix
        gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(M); // workspace for eigenvectors/values
        gsl_vector * eval = gsl_vector_alloc(M);   // eigenvalues
        gsl_matrix * evec = gsl_matrix_alloc(M,M); // eigenvectors
        //
        // gsl_eigen_symmv(pds_covm, eval, evec, w); // NOTE: The diagonal and lower triangular part of A are destroyed during the computation
        gsl_eigen_symmv(&pds_covm_view.matrix, eval, evec, w); // NOTE: The diagonal and lower triangular part of A are destroyed during the computation
        //
        double sigma[M];
        for (size_t i=0; i<M; ++i) {
            // ORSA_DEBUG("eval[%i] = %g",i,gsl_vector_get(eval,i));
            if (gsl_vector_get(eval,i) == 0.0) {
                ORSA_ERROR("problems with the covariance matrix: null eigenvalue found.");
            }
            sigma[i] = sqrt(fabs(gsl_vector_get(eval,i)));
            // ORSA_DEBUG("sigma[%i] = %g",i,sigma[i]);
        }
        gsl_vector * sampleCoeff_x  = gsl_vector_alloc(M);
        gsl_vector * sampleCoeff_y  = gsl_vector_alloc(M); 
        
        std::vector< std::vector<mpf_class> > uniformShape_norm_C;
        std::vector< std::vector<mpf_class> > uniformShape_norm_S;
        {
            uniformShape_norm_C.resize(SH_degree+1);
            uniformShape_norm_S.resize(SH_degree+1);
            for (size_t l=0; l<=SH_degree; ++l) {
                uniformShape_norm_C[l].resize(l+1);
                uniformShape_norm_S[l].resize(l+1);
                for (size_t m=0; m<=l; ++m) {
                    uniformShape_norm_C[l][m] = 0.0;
                    uniformShape_norm_S[l][m] = 0.0;
                }
            }
            
            if (massDistribution.get() != 0) {
                CubicChebyshevMassDistribution::CoefficientType md_uS_coeff;
                CubicChebyshevMassDistribution::resize(md_uS_coeff,0);
                md_uS_coeff[0][0][0] = 1.0;
                osg::ref_ptr<CubicChebyshevMassDistribution> md_uS =
                    new CubicChebyshevMassDistribution(md_uS_coeff,
                                                       1.0, // densityScale = bulk density
                                                       plateModelR0,
                                                       0); // massDistribution->layerData.get());
                orsa::Cache<orsa::Vector> CM = sampled_CM;
                CM.lock();
                CCMD2SH(CM,
                        uniformShape_norm_C,
                        uniformShape_norm_S,
                        SH_degree, // gravityData->degree,
                        si.get(),
                        md_uS,
                        plateModelR0,
                        gravityData->R0);
                
                // const double layerMassFraction = massDistribution->layerData->totalExcessMass() / (GM/orsa::Unit::G());
                /* for (size_t l=0; l<=SH_degree; ++l) {
                   for (size_t m=0; m<=l; ++m) {
                   uniformShape_norm_C[l][m] *= layerMassFraction;
                   uniformShape_norm_S[l][m] *= layerMassFraction;
                   }
                   }
                */
                
            }
        }
        for (size_t l=0; l<=SH_degree; ++l) {
            for (size_t m=0; m<=l; ++m) {
                ORSA_DEBUG("uniformShape_norm_C[%i][%i] = %Fg",l,m,uniformShape_norm_C[l][m].get_mpf_t());
                if (m != 0) ORSA_DEBUG("uniformShape_norm_S[%i][%i] = %Fg",l,m,uniformShape_norm_S[l][m].get_mpf_t());
            }
        }
        
        std::vector< std::vector<mpf_class> > layerData_norm_C;
        std::vector< std::vector<mpf_class> > layerData_norm_S;
        {
            layerData_norm_C.resize(SH_degree+1);
            layerData_norm_S.resize(SH_degree+1);
            for (size_t l=0; l<=SH_degree; ++l) {
                layerData_norm_C[l].resize(l+1);
                layerData_norm_S[l].resize(l+1);
                for (size_t m=0; m<=l; ++m) {
                    layerData_norm_C[l][m] = 0.0;
                    layerData_norm_S[l][m] = 0.0;
                }
            }
            
            if (massDistribution.get() != 0) {
                if (massDistribution->layerData.get() != 0) {
                    CubicChebyshevMassDistribution::CoefficientType md_lD_coeff;
                    CubicChebyshevMassDistribution::resize(md_lD_coeff,0);
                    md_lD_coeff[0][0][0] = 0;
                    osg::ref_ptr<CubicChebyshevMassDistribution> md_lD =
                        new CubicChebyshevMassDistribution(md_lD_coeff,
                                                           0.0,    
                                                           plateModelR0,
                                                           massDistribution->layerData.get());
                    orsa::Cache<orsa::Vector> CM = sampled_CM;
                    CM.lock();
                    CCMD2SH(CM,
                            layerData_norm_C,
                            layerData_norm_S,
                            SH_degree, // gravityData->degree,
                            si.get(),
                            md_lD,
                            plateModelR0,
                            gravityData->R0);
                    
                    const double layerMassFraction = massDistribution->layerData->totalExcessMass() / (GM/orsa::Unit::G());
                    for (size_t l=0; l<=SH_degree; ++l) {
                        for (size_t m=0; m<=l; ++m) {
                            layerData_norm_C[l][m] *= layerMassFraction;
                            layerData_norm_S[l][m] *= layerMassFraction;
                        }
                    }
                    
                }
            }
        }

        for (size_t i=0; i<M; ++i) {
            gsl_vector_set(sampleCoeff_x,i,orsa::GlobalRNG::instance()->rng()->gsl_ran_gaussian(sigma[i]));
        }
        //
        gsl_blas_dgemv(CblasNoTrans,1.0,evec,sampleCoeff_x,0.0,sampleCoeff_y);
        //
        for (size_t i=0; i<M; ++i) {
            
            // correction due to layers
            orsa::Cache<double> layer_coeff;
            if (massDistribution.get() == 0) {
                layer_coeff = 0.0;
            } else if (massDistribution->layerData.get() == 0) {
                layer_coeff = 0.0;
            } else {
                const QString ref_key = mod_gravityData_key(gravityData.get(),i);
                if (ref_key == "GM") {
                    // ORSA_DEBUG("found: [%s]",ref_key.toStdString().c_str());
                    layer_coeff = massDistribution->layerData->totalExcessMass()*orsa::Unit::G();
                } else {
                    for (size_t l=1; l<=SH_degree; ++l) {
                        for (size_t m=0; m<=l; ++m) {
                            if (orsaPDS::RadioScienceGravityData::keyC(l,m) == ref_key) {
                                // ORSA_DEBUG("found: [%s] for l=%i, m=%i",ref_key.toStdString().c_str(),l,m);
                                layer_coeff = layerData_norm_C[l][m].get_d();
                            } else if (orsaPDS::RadioScienceGravityData::keyS(l,m) == ref_key) {
                                // ORSA_DEBUG("found: [%s] for l=%i, m=%i",ref_key.toStdString().c_str(),l,m);
                                layer_coeff = layerData_norm_S[l][m].get_d();
                            }                       
                        }
                    }
                }
            }
            if (!layer_coeff.isSet()) {
                ORSA_DEBUG("problems...");
                exit(0);
            }
            // choose here if sampling or using nominal value
            // nominal (no covariance sampling)
            const double sampled_coeff = gsl_vector_get(pds_coeff,i) - layer_coeff;
            // use covariance sampling
            // const double sampled_coeff = gsl_vector_get(pds_coeff,i) - layer_coeff + gsl_vector_get(sampleCoeff_y,i);
            
            
            // gsl_vector_set(sampleCoeff_y,i,sampled_coeff);
            // gsl_vector_set(sh,z_sh,gsl_vector_get(sampleCoeff_y,z_sh));
            //
            gsl_vector_set(sh,i,sampled_coeff);
            
            // get again pds_covm because old one has been destroyed by the call to gsl_eigen_symmv
            gsl_matrix * pds_covm  = mod_gravityData_getCovarianceMatrix(gravityData.get());
            
            ORSA_DEBUG("%7s = %12.6g [sampled] = %12.6g [layers] + %12.6g   nominal: %+12.6g   delta: %+12.6g   sigma: %12.6g",
                       mod_gravityData_key(gravityData.get(),i).toStdString().c_str(),
                       gsl_vector_get(sh,i)+(*layer_coeff),
                       (*layer_coeff),
                       gsl_vector_get(sh,i),
                       mod_gravityData_getCoeff(gravityData.get(),mod_gravityData_key(gravityData.get(),i)),
                       gsl_vector_get(sh,i)+(*layer_coeff)-mod_gravityData_getCoeff(gravityData.get(),mod_gravityData_key(gravityData.get(),i)),
                       sqrt(gsl_matrix_get(pds_covm,i,i)));
            
            gsl_matrix_free(pds_covm);
        }
        
        // solving here!
        gsl_blas_dgemv(CblasNoTrans,1.0,pseudoInvA,sh,0.0,cT);
        
        gsl_vector * cT0 = gsl_vector_calloc(N);
        gsl_vector_memcpy(cT0,cT);
        
        if (0) {
            
            // trying simulated annealing approach
            
            gsl_rng * rng = ::gsl_rng_alloc(gsl_rng_gfsr4);
            const int randomSeed = time(NULL)*getpid();
            ::gsl_rng_set(rng,randomSeed);
            ORSA_DEBUG("simulated annealing random seed: %d",randomSeed);
            
            std::vector<orsa::Vector> rv;
            //
            {
                const bool storeSamplePoints = false; // saving the points in rv
                osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape =
                    new orsa::RandomPointsInShape(shapeModel,
                                                  0,
                                                  numSamplePoints,
                                                  storeSamplePoints);
                orsa::Vector v;
                randomPointsInShape->reset();
                while (randomPointsInShape->get(v)) {
                    rv.push_back(v);
                }
            }
            
            SIMAN_xp x0;
            x0.R0_plate   = plateModelR0;
            x0.R0_gravity = gravityData->R0;
            x0.bulkDensity = bulkDensity;
            x0.rv = rv;
            x0.SH_degree = SH_degree;
            x0.T_degree = T_degree;
            x0.T_size = T_size;
            x0.cT0 = cT0;
            x0.uK = &uK[0];
            x0.uK_size = N-M;
            x0.factor.resize(x0.uK_size);
            x0.minimumDensity = orsa::FromUnits(orsa::FromUnits(2.50,orsa::Unit::GRAM),orsa::Unit::CM,-3);
            x0.maximumDensity = orsa::FromUnits(orsa::FromUnits(8.00,orsa::Unit::GRAM),orsa::Unit::CM,-3);
            x0.penaltyThreshold = 1.00;
            
            if (massDistribution.get() != 0) {
                x0.layerData = massDistribution->layerData;
            } else {
                x0.layerData = 0;
            }
            
            // fix value of x0.factor[]  
            
            if (have_CCMDF_file) {
                
                // using the input CCMDF
                const size_t cT_CCMDF_degree = massDistribution->coeff.size()-1;
                const size_t cT_CCMDF_size   = CubicChebyshevMassDistribution::totalSize(cT_CCMDF_degree);
                size_t Tx,Ty,Tz;
                for (size_t b=0; b<x0.uK_size; ++b) {
                    x0.factor[b] = 0.0;
                    for (size_t s=0; s<N; ++s) {
                        double cT_CCMDF = 0.0; // default value
                        if (s < cT_CCMDF_size) {
                            CubicChebyshevMassDistribution::triIndex(Tx,Ty,Tz,s);
                            cT_CCMDF = massDistribution->coeff[Tx][Ty][Tz];
                        }
                        x0.factor[b] += (cT_CCMDF-gsl_vector_get(cT0,s))*gsl_vector_get(uK[b],s);
                    }
                    ORSA_DEBUG("factor[%03i] = %g",b,x0.factor[b]);
                }
                
            } else {
                // get as close as possible to cT = {1,0,0,0,0...} = constant density
                // project (1,0,0,0..) - cT0 along uK_b
                for (size_t b=0; b<x0.uK_size; ++b) {
                    x0.factor[b] = 0.0;
                    for (size_t s=0; s<N; ++s) {
                        if (s==0) {
                            // first element of target cT = 1
                            x0.factor[b] += (1.0-gsl_vector_get(cT0,s))*gsl_vector_get(uK[b],s);
                        } else {
                            // all other elements of target cT = 0
                            x0.factor[b] += (0.0-gsl_vector_get(cT0,s))*gsl_vector_get(uK[b],s);
                        }
                    }
                    ORSA_DEBUG("factor[%03i] = %g",b,x0.factor[b]);
                }
            }
            
            gsl_siman_solve(rng, &x0, E1, S1, M1, P1,
                            SIMAN_copy, SIMAN_copy_construct, SIMAN_destroy,
                            0, params);
            
        }
        
        if (1) {
            
            // moving layers...
            
            const size_t uK_size = N-M;
            std::vector<double> factor;
            factor.resize(uK_size);
            const double minimumDensity = orsa::FromUnits(orsa::FromUnits(2.50,orsa::Unit::GRAM),orsa::Unit::CM,-3);
            const double maximumDensity = orsa::FromUnits(orsa::FromUnits(8.00,orsa::Unit::GRAM),orsa::Unit::CM,-3);
            const double penaltyThreshold = 1.00;
            
            std::vector<orsa::Vector> rv;
            //
            {
                const bool storeSamplePoints = false; // saving the points in rv
                osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape =
                    new orsa::RandomPointsInShape(shapeModel,
                                                  0,
                                                  numSamplePoints,
                                                  storeSamplePoints);
                orsa::Vector v;
                randomPointsInShape->reset();
                while (randomPointsInShape->get(v)) {
                    rv.push_back(v);
                }
            }
            
            // fix value of x0.factor[]  
            
            /* 
               if (have_CCMDF_file) {
               // using the input CCMDF
               const size_t cT_CCMDF_degree = massDistribution->coeff.size()-1;
               const size_t cT_CCMDF_size   = CubicChebyshevMassDistribution::totalSize(cT_CCMDF_degree);
               size_t Tx,Ty,Tz;
               for (size_t b=0; b<uK_size; ++b) {
               factor[b] = 0.0;
               for (size_t s=0; s<N; ++s) {
               double cT_CCMDF = 0.0; // default value
               if (s < cT_CCMDF_size) {
               CubicChebyshevMassDistribution::triIndex(Tx,Ty,Tz,s);
               cT_CCMDF = massDistribution->coeff[Tx][Ty][Tz];
               }
               factor[b] += (cT_CCMDF-gsl_vector_get(cT0,s))*gsl_vector_get(uK[b],s);
               }
               ORSA_DEBUG("factor[%03i] = %g",b,factor[b]);
               }
               } else 
            */
            {
                //
                // get as close as possible to cT = {1,0,0,0,0...} = constant density
                // project (1,0,0,0..) - cT0 along uK_b
                for (size_t b=0; b<uK_size; ++b) {
                    factor[b] = 0.0;
                    for (size_t s=0; s<N; ++s) {
                        if (s==0) {
                            // first element of target cT = 1
                            factor[b] += (1.0-gsl_vector_get(cT0,s))*gsl_vector_get(uK[b],s);
                        } else {
                            // all other elements of target cT = 0
                            factor[b] += (0.0-gsl_vector_get(cT0,s))*gsl_vector_get(uK[b],s);
                        }
                    }
                    ORSA_DEBUG("factor[%03i] = %g",b,factor[b]);
                }
            }
            
            gsl_vector * cT = gsl_vector_alloc(T_size);
            gsl_vector_memcpy(cT,cT0);
            for (size_t b=0; b<uK_size; ++b) {
                for (size_t j=0; j<T_size; ++j) {
                    gsl_vector_set(cT,j,gsl_vector_get(cT,j)+factor[b]*gsl_vector_get(uK[b],j));
                }
            }
            
            CubicChebyshevMassDistribution::CoefficientType coeff;
            CubicChebyshevMassDistribution::resize(coeff,T_degree); 
            
            for (unsigned int i=0; i<=T_degree; ++i) {
                for (unsigned int j=0; j<=T_degree-i; ++j) {
                    for (unsigned int k=0; k<=T_degree-i-j; ++k) {
                        if (i+j+k<=T_degree) {
                            const size_t index = CubicChebyshevMassDistribution::index(i,j,k);
                            coeff[i][j][k] = gsl_vector_get(cT,index);
                        }
                    }            
                }
            }
            
            osg::ref_ptr<CubicChebyshevMassDistribution> newMassDistribution =
                new CubicChebyshevMassDistribution(coeff,bulkDensity,plateModelR0,massDistribution->layerData);
            
            osg::ref_ptr< orsa::Statistic<double> > stat = new orsa::Statistic<double>;
            orsa::Cache<double> minDensity, maxDensity;
            std::vector<double> dv;
            dv.resize(rv.size());
            for (size_t k=0; k<rv.size(); ++k) {
                dv[k] = newMassDistribution->density(rv[k]);
                stat->insert(dv[k]);
                minDensity.setIfSmaller(dv[k]);
                maxDensity.setIfLarger(dv[k]);
            }
            
            const double averageDensity = stat->average();
            
            const double penalty =
                MassDistributionPenalty(rv,
                                            dv,
                                        newMassDistribution.get());
            
#warning find a better way to compute average density? (based on simplexIntegral and not on randomPointsInShape)
            
            ORSA_DEBUG("[density] min: %+6.2f max: %+6.2f avg: %+6.2f [g/cm^3]   penalty: %g",
                       orsa::FromUnits(orsa::FromUnits(minDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
                       orsa::FromUnits(orsa::FromUnits(maxDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
                       orsa::FromUnits(orsa::FromUnits(averageDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
                       penalty);
            
            if ( (minDensity >= minimumDensity) &&
                 (maxDensity <= maximumDensity) &&
                 (penalty <= penaltyThreshold) ) {
                // another quick output...
#warning pass filename as parameter...
                CubicChebyshevMassDistributionFile::CCMDF_data data;
                data.minDensity = minDensity;
                data.maxDensity = maxDensity;
                data.deltaDensity = maxDensity-minDensity;
                data.penalty = penalty;
                data.densityScale = bulkDensity;
                data.R0 = plateModelR0;
                data.SH_degree = SH_degree;
                data.coeff = coeff;
                data.layerData = newMassDistribution->layerData;
                CubicChebyshevMassDistributionFile::append(data,"CCMDF.out");
            }
            
        }
        
    }
    
    // free GSL stuff
    gsl_vector_free(pds_coeff);
    gsl_matrix_free(pds_covm);
    gsl_matrix_free(pds_inv_covm);
    gsl_matrix_free(cT2sh);
    
#warning call gsl_*_free as needed...
    
    return 0;
}

