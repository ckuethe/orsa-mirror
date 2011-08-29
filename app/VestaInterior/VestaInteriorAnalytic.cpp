#include "VestaInteriorAnalytic.h"

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

/*** CHOOSE ONE ***/
// typedef double simplex_T;
// typedef mpf_class simplex_T;
typedef dd_real simplex_T;
// typedef qd_real simplex_T;

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;


int main(int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    // QD
    unsigned int oldcw;
    fpu_fix_start(&oldcw);
    
    //ORSA_DEBUG("current mpf precision: %i",mpf_get_default_prec());
    mpf_set_default_prec(128);
    // mpf_set_default_prec(256);
    // mpf_set_default_prec(512);
    // ORSA_DEBUG("updated mpf precision: %i",mpf_get_default_prec());
    
    if (argc != 13) {
        printf("Usage: %s <RadioScienceGravityFile> <plate-model-file> <plate-model-R0_km> <gravity-degree> <polynomial-degree> <CM-x_km> <CM-y_km> <CM-z_km> <CM-sigma-x_km> <CM-sigma-y_km> <CM-sigma-z_km> <num-sample-points>\n",argv[0]);
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
    
    // safer over NFS
    sqlite3_vfs_register(sqlite3_vfs_find("unix-dotfile"), 1);
    
    const std::string SQLiteDBFileName = getSqliteDBFileName(plateModelFile,plateModelR0);
    
    /* // test
       {
       for (size_t n=0; n<=10; ++n) {
       const std::vector<mpz_class> & coeff = orsa::ChebyshevTcoeff(n);
       gmp_printf("cT[n=%03i]: ",n);
       for (size_t k=0; k<=n; ++k) {
       gmp_printf("%+9Zi ",coeff[k].get_mpz_t());
       }
       gmp_printf("\n");
       }
       }
    */
    /* {
       for (size_t n=0; n<=10; ++n) {
       const std::vector<double> & ext = orsa::ChebyshevTextrema(n);
       gmp_printf("cText[n=%03i]: ",n);
       for (size_t k=0; k<=n; ++k) {
       gmp_printf("%+12.9f ",ext[k]);
       }
       gmp_printf("\n");
       }
       }
    */
    /* {
       orsa::triIndex_mpq tri_integral;
       orsa::triIndex_d   tri_plain;
       orsa::triIndex_d   tri_norm;
       for (size_t l=0; l<=10; ++l) {
       for(size_t m=0; m<=l; ++m) {
       {
       // C
       tri_integral = orsa::conversionCoefficients_C_integral(l,m);
       tri_plain    = orsa::conversionCoefficients_C_plain(l,m);
       tri_norm     = orsa::conversionCoefficients_C_norm(l,m);
       for (size_t ti=0; ti<=l; ++ti) {
       for (size_t tj=0; tj<=l-ti; ++tj) {
       for (size_t tk=0; tk<=l-ti-tj; ++tk) {
       if (tri_integral[ti][tj][tk] != 0) {
       if (tri_integral[ti][tj][tk].get_den() == 1) {
       ORSA_DEBUG("C[%i][%i] += %Zi N[%i][%i][%i]",l,m,
       tri_integral[ti][tj][tk].get_num().get_mpz_t(),
       // tri_integral[ti][tj][tk].get_den().get_mpz_t(),
       ti,tj,tk);
       } else {
       ORSA_DEBUG("C[%i][%i] += %Zi/%Zi N[%i][%i][%i]",l,m,
       tri_integral[ti][tj][tk].get_num().get_mpz_t(),
       tri_integral[ti][tj][tk].get_den().get_mpz_t(),
       ti,tj,tk);
       }
       ORSA_DEBUG("C[%i][%i] += %g N[%i][%i][%i]",l,m,
       tri_plain[ti][tj][tk],
       ti,tj,tk);
       ORSA_DEBUG("C[%i][%i] += %g N[%i][%i][%i]",l,m,
       tri_norm[ti][tj][tk],
       ti,tj,tk);
       }
       }
       }
       }
       }
       {
       // S
       tri_integral = orsa::conversionCoefficients_S_integral(l,m);
       tri_plain    = orsa::conversionCoefficients_S_plain(l,m);
       tri_norm     = orsa::conversionCoefficients_S_norm(l,m);
       for (size_t ti=0; ti<=l; ++ti) {
       for (size_t tj=0; tj<=l-ti; ++tj) {
       for (size_t tk=0; tk<=l-ti-tj; ++tk) {
       if (tri_integral[ti][tj][tk] != 0) {
       if (tri_integral[ti][tj][tk].get_den() == 1) {
       ORSA_DEBUG("S[%i][%i] += %Zi N[%i][%i][%i]",l,m,
       tri_integral[ti][tj][tk].get_num().get_mpz_t(),
       // tri_integral[ti][tj][tk].get_den().get_mpz_t(),
       ti,tj,tk);
       } else {
       ORSA_DEBUG("S[%i][%i] += %Zi/%Zi N[%i][%i][%i]",l,m,
       tri_integral[ti][tj][tk].get_num().get_mpz_t(),
       tri_integral[ti][tj][tk].get_den().get_mpz_t(),
       ti,tj,tk);
       }
       ORSA_DEBUG("S[%i][%i] += %g N[%i][%i][%i]",l,m,
       tri_plain[ti][tj][tk],
       ti,tj,tk);
       ORSA_DEBUG("S[%i][%i] += %g N[%i][%i][%i]",l,m,
       tri_norm[ti][tj][tk],
       ti,tj,tk);
       }
       }
       }
       }
       }
       }
       }
       }
    */
    
    // test specific cases, for debug purposes only!
    // orsa::GlobalRNG::randomSeed = -800402816;
    
    osg::ref_ptr<orsaPDS::RadioScienceGravityData> gravityData = new orsaPDS::RadioScienceGravityData;
    orsaPDS::RadioScienceGravityFile::read(gravityData.get(),radioScienceGravityFile,512,1518);
    
    if (0) {
        
        ORSA_DEBUG("testing writing of gravity data file...");
        
        // test
        std::string testGravityFile = "gra.dat";
        orsaPDS::RadioScienceGravityFile::write(gravityData.get(),testGravityFile,512,1518);
        
        osg::ref_ptr<orsaPDS::RadioScienceGravityData> testGravityData = new orsaPDS::RadioScienceGravityData;
        orsaPDS::RadioScienceGravityFile::read(testGravityData.get(),testGravityFile,512,1518);
        
        {
            gsl_vector *  pds_coeff = gravityData->getCoefficientVector();
            gsl_vector * test_coeff = testGravityData->getCoefficientVector();
            for (size_t i=0; i<gravityData->numberOfCoefficients; ++i) {
                const double delta =
                    gsl_vector_get(test_coeff,i) -
                    gsl_vector_get( pds_coeff,i);
                const double perc = 100*fabs(delta / gsl_vector_get(pds_coeff,i));
                if (delta != 0.0) {
                    ORSA_DEBUG("i: %i %+12.6e %9.3e %%",i,delta,perc);
                }
            }
        }
        
        {
            gsl_matrix *  pds_covm = gravityData->getCovarianceMatrix();
            gsl_matrix * test_covm = testGravityData->getCovarianceMatrix();
            for (size_t i=0; i<gravityData->numberOfCoefficients; ++i) {
                for (size_t j=0; j<=i; ++j) {
                    const double delta =
                        gsl_matrix_get(test_covm,i,j) -
                        gsl_matrix_get( pds_covm,i,j);
                    const double perc = 100*fabs(delta / gsl_matrix_get(pds_covm,i,j));
                    if (delta != 0.0) {
                        ORSA_DEBUG("i: %i j: %i %+12.6e %9.3e %%",i,j,delta,perc);
                    }
                }
            }
        }
        
        ORSA_DEBUG("testing writing of gravity data file... done.");
    }
    
    /* osg::ref_ptr<orsaPDS::RadioScienceGravityFile> pds =
    // new orsaPDS::RadioScienceGravityFile("JGDAWN20SIMA.DAT",512,1518);
    new orsaPDS::RadioScienceGravityFile(radioScienceGravityFile,512,1518);
    */
    
    gsl_vector * pds_coeff    = gravityData->getCoefficientVector();
    gsl_matrix * pds_covm     = gravityData->getCovarianceMatrix();
    gsl_matrix * pds_inv_covm = gravityData->getInverseCovarianceMatrix();
    
    /* if (0) {
       gsl_matrix * identity = gsl_matrix_alloc(gravityData->numberOfCoefficients,gravityData->numberOfCoefficients);
       // gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,pds_covm,pds_inv_covm,0.0,identity);
       gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,pds_inv_covm,pds_covm,0.0,identity);
       for (size_t i=0; i<gravityData->numberOfCoefficients; ++i) {
       for (size_t j=0; j<gravityData->numberOfCoefficients; ++j) {
       ORSA_DEBUG("identity[%03i][%03i] = %+20.9f",i,j,gsl_matrix_get(identity,i,j));
       }
       }
       }
    */
    
    /* osg::ref_ptr<VestaShape> shape = new VestaShape;
       if (!shape->read("vesta_thomas.dat")) {
       ORSA_ERROR("problems encountered while reading shape file...");
       exit(0);
       }
    */
    
    osg::ref_ptr<GaskellPlateModel> shapeModel = new GaskellPlateModel;
    if (!shapeModel->read(plateModelFile)) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    osg::ref_ptr<SimplexIntegration<simplex_T> > si = new SimplexIntegration<simplex_T>(shapeModel.get(), plateModelR0, SQLiteDBFileName);
    si->reserve(polynomialDegree);
    
    const size_t SH_degree = gravityDegree; // shperical harmonics degree
    const size_t  T_degree = polynomialDegree; // chebyshev polinomials degree
    
#warning which GM value to use? gravityData->GM  OR gravityData->getCoeff("GM") ??
    const double GM = gravityData->GM; 
    
    // sh size: (l+1)^2 +1; +1 due to GM factor; -4 because C00, C10, C11, S11 are missing
    const size_t  SH_size = (SH_degree+1)*(SH_degree+1)+1-4;
    const size_t ijk_size = CubicChebyshevMassDistribution::totalSize(SH_degree);
    const size_t   T_size = CubicChebyshevMassDistribution::totalSize( T_degree);
    
    ORSA_DEBUG("SH_size: %d  T_size: %d",SH_size,T_size);
    
    if (T_size <= SH_size) {
        ORSA_DEBUG("this method works only when the problem is under-determined, exiting");
        exit(0);
    }
    
    gsl_matrix *      sh2ijk = gsl_matrix_calloc(SH_size,ijk_size);
    //
    gsl_matrix *   pq_sh2ijk =  gsl_matrix_alloc(SH_size,ijk_size);
    gsl_matrix *   nu_sh2ijk =  gsl_matrix_alloc(SH_size,ijk_size);
    //
    gsl_matrix * identity_sh = gsl_matrix_calloc(SH_size,SH_size);
    for (size_t z_sh=0; z_sh<SH_size; ++z_sh) {
        gsl_matrix_set(identity_sh,z_sh,z_sh,1.0);
    }

    
    // C_lm coefficients
    //
    for (int l=0; l<=(int)SH_degree; ++l) {

#warning remember: skipping all l==1 terms
        if (l==1) continue;
        
        for (int m=0; m<=l; ++m) {
            
            // double pq_factor=0;
            // double pq_factor_uncertainty=0;
            //
            // reset pq_sh2ijk to zero
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, 0.0, pq_sh2ijk);
            
            // integer division in limits
            for (int p=0;p<=(l/2);++p) {
                for (int q=0;q<=(m/2);++q) {
	  
                    // double nu_factor=0;
                    // double nu_factor_uncertainty=0;
                    //
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, nu_sh2ijk, 0.0, nu_sh2ijk);
                    
                    for (int nu_x=0; nu_x<=p; ++nu_x) {
                        for (int nu_y=0; nu_y<=(p-nu_x); ++nu_y) {
                            
                            const int M_i = m-2*q+2*nu_x;
                            const int M_j = 2*q+2*nu_y;
                            const int M_k = l-m-2*nu_x-2*nu_y;
                            
                            if (M_i+M_j+M_k!=l) {
                                ORSA_DEBUG("WARNING!!! ***********************");
                            }
	      
                            if ( (M_i>=0) && 
                                 (M_j>=0) && 
                                 (M_k>=0) && 
                                 (M_i+M_j+M_k==l) ) {
		
                                // ORSA_DEBUG("requesting M[%i][%i][%i]...   l: %i",M_i, M_j, M_k, l);
		
                                const double nu_factor_base =
                                    (orsa::factorial(p).get_d() / (orsa::factorial(nu_x).get_d()*orsa::factorial(nu_y).get_d()*orsa::factorial(p-nu_x-nu_y).get_d()));
		
                                // nu_factor += nu_factor_base * pm->M(M_i,M_j,M_k);
                                // nu_factor_uncertainty += nu_factor_base * pm->M_uncertainty(M_i,M_j,M_k);

                                // z_sh is 0 if l==0, it's the GM term
                                const size_t z_sh  = (l==0 ? 0 :  gravityData->index(orsaPDS::RadioScienceGravityData::keyC(l,m)));
                                const size_t z_ijk = CubicChebyshevMassDistribution::index(M_i,M_j,M_k);
                                // gsl_matrix_set(nu_sh2ijk,z_sh,z_ijk,gsl_matrix_get(nu_sh2ijk,z_sh,z_ijk)+nu_factor_base);
                                gsl_matrix_set(nu_sh2ijk,z_sh,z_ijk,nu_factor_base);
                                
                            }
                        }
                    }
                    // need this because using M(i,j,k) instead of N(i,j,k)
                    // nu_factor /= int_pow(R0,l);
                    // nu_factor_uncertainty /= int_pow(R0,l);

#warning RESTORE?
                    // gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, nu_sh2ijk, 1.0/int_pow(R0,l), nu_sh2ijk);
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, nu_sh2ijk, 1.0, nu_sh2ijk);

                    /* 
                       const double pq_factor_base = 
                       orsa::power_sign(p+q) *
                       orsa::binomial(l,p).get_d() *
                       orsa::binomial(2*l-2*p,l).get_d() *
                       orsa::binomial(m,2*q).get_d() *
                       orsa::pochhammer(mpz_class(l-m-2*p+1),m).get_d();
                    */
                    const double pq_factor_base = 
                        mpz_class(orsa::power_sign(p+q) *
                                  orsa::binomial(l,p) *
                                  orsa::binomial(2*l-2*p,l) *
                                  orsa::binomial(m,2*q) *
                                  orsa::pochhammer(mpz_class(l-m-2*p+1),m)).get_d();
                    
                    /* pq_factor += 
                       pq_factor_base *
                       nu_factor;
                       
                       pq_factor_uncertainty += 
                       pq_factor_base *
                       nu_factor_uncertainty;
                    */
                    
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, nu_sh2ijk, pq_factor_base, nu_sh2ijk);
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, identity_sh, nu_sh2ijk, 1.0, pq_sh2ijk);
                }
            }
      
            // pq_factor /= int_pow(2.0,l);
            // pq_factor_uncertainty /= int_pow(2.0,l);
            //
            
            // #warning RESTORE
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, 1.0/orsa::int_pow(2.0,l), pq_sh2ijk);
            // gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, 1.0, pq_sh2ijk);
            
            // ORSA_DEBUG("PaulMoment::normalization(%i,%i) = %g",l,m,PaulMoment::normalization(l,m));
            // normalization
            // gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, PaulMoment::normalization(l,m), pq_sh2ijk);
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, orsa::normalization_integralToNormalizedSphericalHarmonics(l,m).get_d(), pq_sh2ijk);


#warning scaling for l=0,m=0 term that is GM in the data...
            if (l==0 && m==0) {
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, GM, pq_sh2ijk);
            }
            
            /* const double C_lm = pq_factor;
               const double C_lm_uncertainty = fabs(pq_factor_uncertainty);
               //      
               const double norm_C_lm = C_lm*PaulMoment::normalization(l,m);
               const double norm_C_lm_uncertainty = fabs(C_lm_uncertainty*PaulMoment::normalization(l,m));
            */
            //
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, identity_sh, pq_sh2ijk, 1.0, sh2ijk);
            
            /* if (verbose) {
               ORSA_DEBUG("     C[%i][%i] = %+16.12f +/- %16.12f",
               l,m,     C_lm, C_lm_uncertainty);
               ORSA_DEBUG("norm_C[%i][%i] = %+16.12f +/- %16.12f   norm: %f",
               l,m,norm_C_lm,norm_C_lm_uncertainty,PaulMoment::normalization(l,m));
               }
            */
            
            /* C[l][m]      = C_lm;
               norm_C[l][m] = norm_C_lm;
            */
        }
    } 

    // S_lm coefficients
    //
    for (int l=2; l<=(int)SH_degree; ++l) {
        for (int m=1; m<=l; ++m) {
            
            // double pq_factor=0;
            // double pq_factor_uncertainty=0;
            //
            // reset pq_sh2ijk to zero
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, 0.0, pq_sh2ijk);
            
            // integer division in limits
            for (int p=0;p<=(l/2);++p) {
                for (int q=0;q<=((m-1)/2);++q) {
	  
                    // double nu_factor=0;
                    // double nu_factor_uncertainty=0;
                    //
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, nu_sh2ijk, 0.0, nu_sh2ijk);
                    
                    for (int nu_x=0; nu_x<=p; ++nu_x) {
                        for (int nu_y=0; nu_y<=(p-nu_x); ++nu_y) {
	      
                            const int M_i = m-2*q-1+2*nu_x;
                            const int M_j = 2*q+1+2*nu_y;
                            const int M_k = l-m-2*nu_x-2*nu_y;
	      
                            if (M_i+M_j+M_k!=l) {
                                ORSA_DEBUG("WARNING!!!");
                            }
	      
                            if ( (M_i>=0) && 
                                 (M_j>=0) && 
                                 (M_k>=0) && 
                                 (M_i+M_j+M_k==l) ) {
		
                                // ORSA_DEBUG("requesting M[%i][%i][%i]...   l: %i",M_i, M_j, M_k, l);
		
                                const double nu_factor_base =
                                    (orsa::factorial(p).get_d() /(orsa::factorial(nu_x).get_d()*orsa::factorial(nu_y).get_d()*orsa::factorial(p-nu_x-nu_y).get_d()));
		
                                // nu_factor += nu_factor_base * pm->M(M_i,M_j,M_k);
                                // nu_factor_uncertainty += nu_factor_base * pm->M_uncertainty(M_i,M_j,M_k);
                                
                                const size_t z_sh  = gravityData->index(orsaPDS::RadioScienceGravityData::keyS(l,m));
                                const size_t z_ijk = CubicChebyshevMassDistribution::index(M_i,M_j,M_k);
                                // gsl_matrix_set(nu_sh2ijk,z_sh,z_ijk,gsl_matrix_get(nu_sh2ijk,z_sh,z_ijk)+nu_factor_base);
                                gsl_matrix_set(nu_sh2ijk,z_sh,z_ijk,nu_factor_base);
                                
                            }
                        }
                    }
                    // need this because using M(i,j,k) instead of N(i,j,k)
                    // nu_factor /= int_pow(R0,l);
                    // nu_factor_uncertainty /= int_pow(R0,l);
	  
#warning RESTORE?
                    // gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, nu_sh2ijk, 1.0/int_pow(R0,l), nu_sh2ijk);
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, nu_sh2ijk, 1.0, nu_sh2ijk);
                    
                    
                    /* const double pq_factor_base = 
                       orsa::power_sign(p+q) *
                       orsa::binomial(l,p).get_d() *
                       orsa::binomial(2*l-2*p,l).get_d() *
                       orsa::binomial(m,2*q+1).get_d() *
                       orsa::pochhammer(mpz_class(l-m-2*p+1),m).get_d();
                    */
                    const double pq_factor_base = 
                        mpz_class(orsa::power_sign(p+q) *
                                  orsa::binomial(l,p) *
                                  orsa::binomial(2*l-2*p,l) *
                                  orsa::binomial(m,2*q+1) *
                                  orsa::pochhammer(mpz_class(l-m-2*p+1),m)).get_d();
                    
                    /* pq_factor += 
                       pq_factor_base *
                       nu_factor;
                       
                       pq_factor_uncertainty += 
                       pq_factor_base *
                       nu_factor_uncertainty;
                    */
                    
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, nu_sh2ijk, pq_factor_base, nu_sh2ijk);
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, identity_sh, nu_sh2ijk, 1.0, pq_sh2ijk);
                }
            }
            //
            // pq_factor /= int_pow(2.0,l);
            // pq_factor_uncertainty /= int_pow(2.0,l);
            //
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, 1.0/orsa::int_pow(2.0,l), pq_sh2ijk);
            
            // ORSA_DEBUG("PaulMoment::normalization(%i,%i) = %g",l,m,PaulMoment::normalization(l,m));
            // normalization
            // gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, PaulMoment::normalization(l,m), pq_sh2ijk);
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, orsa::normalization_integralToNormalizedSphericalHarmonics(l,m).get_d(), pq_sh2ijk);
            
            /* const double S_lm = pq_factor;
               const double S_lm_uncertainty = fabs(pq_factor_uncertainty);
               //      
               const double norm_S_lm = S_lm*PaulMoment::normalization(l,m);
               const double norm_S_lm_uncertainty = fabs(S_lm_uncertainty*PaulMoment::normalization(l,m));
            */
            
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, identity_sh, pq_sh2ijk, 1.0, sh2ijk);
            
            /* if (verbose) {
               ORSA_DEBUG("     S[%i][%i] = %+16.12f +/- %16.12f",
               l,m,     S_lm, S_lm_uncertainty);
               ORSA_DEBUG("norm_S[%i][%i] = %+16.12f +/- %16.12f   norm: %f",
               l,m,norm_S_lm,norm_S_lm_uncertainty,PaulMoment::normalization(l,m));
               }
            */
            
            /* S[l][m]      = S_lm;
               norm_S[l][m] = norm_S_lm;
            */
        }
    }
    
#warning check if there is any ROTATION between reference systems
    
#warning FIX GM entry!!
    
#warning re-include factor 1/R0^l ?
    
#warning re-include normalization for Ckm Slm
    
#warning maybe rescale GM entries to 1.0 to have more homogeneous entries?
    
#warning CHECK ALL NORMALIZATIONS!
    
#warning check code for HIGH degree, might have to rewrite linear algebra...
    
    for (size_t z_sh=0; z_sh<SH_size; ++z_sh) {
        for (size_t z_ijk=0; z_ijk<ijk_size; ++z_ijk) {
            if (gsl_matrix_get(sh2ijk,z_sh,z_ijk)!=0.0) {
                size_t nx,ny,nz;
                CubicChebyshevMassDistribution::triIndex(nx,ny,nz,z_ijk);                
                ORSA_DEBUG("sh2ijk[%03i][%03i] = %+10.6f [%s -> N[%02i][%02i][%02i]]",
                           z_sh,z_ijk,gsl_matrix_get(sh2ijk,z_sh,z_ijk),gravityData->key(z_sh).toStdString().c_str(),nx,ny,nz);
            }
        }
    }
    
    const orsa::Vector sampled_CM(CM_x+orsa::GlobalRNG::instance()->rng()->gsl_ran_gaussian(CM_sx),  
                                  CM_y+orsa::GlobalRNG::instance()->rng()->gsl_ran_gaussian(CM_sy),
                                  CM_z+orsa::GlobalRNG::instance()->rng()->gsl_ran_gaussian(CM_sz));
    
    const double volume = si->getIntegral(0,0,0)*orsa::cube(plateModelR0);
    
    const double bulkDensity = GM/orsa::Unit::G()/volume;
    const double bulkDensity_gcm3 = orsa::FromUnits(orsa::FromUnits(bulkDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3);
    
    ORSA_DEBUG("bulkDensity: %g",bulkDensity);
    
#warning fix how bulkDensity is used in this code...
    
    gsl_matrix * ijk2cT = gsl_matrix_calloc(ijk_size,T_size);
    
    const double radiusCorrectionRatio = plateModelR0/gravityData->R0;
    
    for (size_t z_cT=0; z_cT<T_size; ++z_cT) {
        
        size_t Tx,Ty,Tz;
        CubicChebyshevMassDistribution::triIndex(Tx,Ty,Tz,z_cT);
        
        const std::vector<mpz_class> & cTx = orsa::ChebyshevTcoeff(Tx);
        const std::vector<mpz_class> & cTy = orsa::ChebyshevTcoeff(Ty);
        const std::vector<mpz_class> & cTz = orsa::ChebyshevTcoeff(Tz);
        
        for (size_t z_ijk=0; z_ijk<ijk_size; ++z_ijk) {
            
            size_t nx,ny,nz;
            CubicChebyshevMassDistribution::triIndex(nx,ny,nz,z_ijk);
            
            double sum = 0;
            // expansion of Chebyshev
            for (size_t cx=0; cx<=Tx; ++cx) {
                if (cTx[cx] == 0) continue;
                for (size_t cy=0; cy<=Ty; ++cy) {
                    if (cTy[cy] == 0) continue;
                    for (size_t cz=0; cz<=Tz; ++cz) {
                        if (cTz[cz] == 0) continue;
                        
                        // binomial expansion of (vec-vec_CM)
                        for (size_t bx=0; bx<=nx; ++bx) {
                            for (size_t by=0; by<=ny; ++by) {
                                for (size_t bz=0; bz<=nz; ++bz) {
                                    
                                    sum += orsa::power_sign(bx+by+bz) *
                                        mpz_class(orsa::binomial(nx,bx) *
                                                  orsa::binomial(ny,by) *
                                                  orsa::binomial(nz,bz) *
                                                  cTx[cx] *
                                                  cTy[cy] *
                                                  cTz[cz]).get_d() *
                                        orsa::int_pow(sampled_CM.getX()/plateModelR0,bx) *
                                        orsa::int_pow(sampled_CM.getY()/plateModelR0,by) *
                                        orsa::int_pow(sampled_CM.getZ()/plateModelR0,bz) *
                                        si->getIntegral(nx-bx+cx,ny-by+cy,nz-bz+cz);
                                }
                            }
                        }
                    }
                }
            }
            sum *= orsa::int_pow(radiusCorrectionRatio,nx+ny+nz);
            
            // gsl_matrix_set(ijk2cT,z_ijk,z_cT,sum);
            // gsl_matrix_set(ijk2cT,z_ijk,z_cT,sum/(volume*orsa::int_pow(plateModelR0,-3)));
            gsl_matrix_set(ijk2cT,z_ijk,z_cT,sum/si->getIntegral(0,0,0));
            
            if (1) {
                size_t nx,ny,nz;
                CubicChebyshevMassDistribution::triIndex(nx,ny,nz,z_ijk);                
                size_t Tx,Ty,Tz;
                CubicChebyshevMassDistribution::triIndex(Tx,Ty,Tz,z_cT);
                ORSA_DEBUG("ijk2cT[%03i][%03i] = %+9.6f [N[%02i][%02i][%02i] -> cT[%i][%i][%i]]",
                           z_ijk,z_cT,gsl_matrix_get(ijk2cT,z_ijk,z_cT),nx,ny,nz,Tx,Ty,Tz);
            }
        }
    }
    
    
    /* for (size_t z_ijk=0; z_ijk<ijk_size; ++z_ijk) {
       for (size_t z_cT=0; z_cT<T_size; ++z_cT) {
       // if (gsl_matrix_get(ijk2cT,z_ijk,z_cT)!=0.0) {
       {
       size_t nx,ny,nz;
       CubicChebyshevMassDistribution::triIndex(nx,ny,nz,z_ijk);                
       size_t Tx,Ty,Tz;
       CubicChebyshevMassDistribution::triIndex(Tx,Ty,Tz,z_cT);
       ORSA_DEBUG("ijk2cT[%03i][%03i] = %+9.6f [N[%02i][%02i][%02i] -> cT[%i][%i][%i]]",
       z_ijk,z_cT,gsl_matrix_get(ijk2cT,z_ijk,z_cT),nx,ny,nz,Tx,Ty,Tz);
       }
       }
       }
    */
    
    gsl_matrix * sh2cT = gsl_matrix_alloc(SH_size,T_size);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, sh2ijk, ijk2cT, 0.0, sh2cT);
    
    for (size_t z_sh=0; z_sh<SH_size; ++z_sh) {
        for (size_t z_cT=0; z_cT<T_size; ++z_cT) {
            size_t Tx,Ty,Tz;
            CubicChebyshevMassDistribution::triIndex(Tx,Ty,Tz,z_cT);
            ORSA_DEBUG("sh2cT[%03i][%03i] = %+9.6f [%s -> cT[%i][%i][%i]]",
                       z_sh,z_cT,gsl_matrix_get(sh2cT,z_sh,z_cT),gravityData->key(z_sh).toStdString().c_str(),Tx,Ty,Tz);
        }
    }
    
    {
        
        // A x = b
        // A = sh2cT
        // x = A^T (A A^T)^(-1) b
        
        const size_t M = SH_size;
        const size_t N =  T_size;
        
        // A
        gsl_matrix * A = gsl_matrix_alloc(M,N);
        
        for (size_t j=0; j<M; ++j) {
            for (size_t k=0; k<N; ++k) {
                gsl_matrix_set(A,j,k,gsl_matrix_get(sh2cT,j,k));
            }
        }
        
        // A^T
        gsl_matrix * AT = gsl_matrix_alloc(N,M);
        
        for (size_t j=0; j<N; ++j) {
            for (size_t k=0; k<M; ++k) {
                gsl_matrix_set(AT,j,k,gsl_matrix_get(sh2cT,k,j));
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
        
        for (size_t j=0; j<M; ++j) {
            for (size_t k=0; k<M; ++k) {
                ORSA_DEBUG("(A A^T)[%03i][%03i] = %+12.6f",j,k,gsl_matrix_get(A_AT,j,k));
            }
        }
        
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
        gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(gravityData->numberOfCoefficients); // workspace for eigenvectors/values
        gsl_vector * eval = gsl_vector_alloc(gravityData->numberOfCoefficients);    // eigenvalues
        gsl_matrix * evec = gsl_matrix_alloc(gravityData->numberOfCoefficients,gravityData->numberOfCoefficients);  // eigenvectors
        //
        gsl_eigen_symmv(pds_covm, eval, evec, w); // NOTE: The diagonal and lower triangular part of A are destroyed during the computation
        //
        double sigma[gravityData->numberOfCoefficients];
        for (size_t i=0; i<gravityData->numberOfCoefficients; ++i) {
            // ORSA_DEBUG("eval[%i] = %g",i,gsl_vector_get(eval,i));
            if (gsl_vector_get(eval,i) == 0.0) {
                ORSA_ERROR("problems with the covariance matrix: null eigenvalue found.");
            }
            sigma[i] = sqrt(fabs(gsl_vector_get(eval,i)));
            // ORSA_DEBUG("sigma[%i] = %g",i,sigma[i]);
        }
        gsl_vector * sampleCoeff_x  = gsl_vector_alloc(gravityData->numberOfCoefficients);
        gsl_vector * sampleCoeff_y  = gsl_vector_alloc(gravityData->numberOfCoefficients); 
        
        for (size_t gen=0; gen<1000; ++gen) {
            
            for (size_t i=0; i<gravityData->numberOfCoefficients; ++i) {
                gsl_vector_set(sampleCoeff_x,i,orsa::GlobalRNG::instance()->rng()->gsl_ran_gaussian(sigma[i]));
            }
            //
            gsl_blas_dgemv(CblasNoTrans,1.0,evec,sampleCoeff_x,0.0,sampleCoeff_y);
            //
            for (size_t i=0; i<gravityData->numberOfCoefficients; ++i) {
                gsl_vector_set(sampleCoeff_y,i,gsl_vector_get(sampleCoeff_y,i)+gsl_vector_get(pds_coeff,i));
            }
            
            for (size_t z_sh=0; z_sh<SH_size; ++z_sh) {
                // gsl_vector_set(sh,z_sh,gravityData->getCoeff(gravityData->key(z_sh)));
                gsl_vector_set(sh,z_sh,gsl_vector_get(sampleCoeff_y,z_sh));
                ORSA_DEBUG("%7s =  %+12.3g [sampled]",gravityData->key(z_sh).toStdString().c_str(),gsl_vector_get(sh,z_sh));
            }

            // solving here!
            gsl_blas_dgemv(CblasNoTrans,1.0,pseudoInvA,sh,0.0,cT);
            
            gsl_vector * cT0 = gsl_vector_calloc(N);
            gsl_vector_memcpy(cT0,cT);


            if (1) {
                
                // trying simulated annealing approach
                
                gsl_rng * rng = ::gsl_rng_alloc(gsl_rng_gfsr4);
                const int randomSeed = time(NULL)*getpid();
                ::gsl_rng_set(rng,randomSeed);
                ORSA_DEBUG("simulated annealing random seed: %d",randomSeed);
                
                // const size_t numSamplePoints = 100000;
                const bool storeSamplePoints = true;
                osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape =
                    new orsa::RandomPointsInShape(shapeModel,
                                                  0,
                                                  numSamplePoints,
                                                  storeSamplePoints);
                
                SIMAN_xp x0;
#warning which R0 to use (both?)
                x0.R0_plate   = plateModelR0;
                x0.R0_gravity = gravityData->R0;
                x0.bulkDensity_gcm3 = bulkDensity_gcm3;
                x0.randomPointsInShape = randomPointsInShape;
                x0.T_degree = T_degree;
                x0.T_size = T_size;
                x0.cT0 = cT0;
                x0.uK = &uK[0];
                x0.uK_size = N-M;
                x0.factor.resize(x0.uK_size);
                for (size_t b=0; b<x0.uK_size; ++b) {
                    x0.factor[b] = 0.0;
                }
                
                gsl_siman_solve(rng, &x0, E1, S1, M1, P1,
                                SIMAN_copy, SIMAN_copy_construct, SIMAN_destroy,
                                0, params);
                
            }
            
        }
    }
    
    // free GSL stuff
    gsl_vector_free(pds_coeff);
    gsl_matrix_free(pds_covm);
    gsl_matrix_free(pds_inv_covm);
    gsl_matrix_free(sh2ijk);
    gsl_matrix_free(pq_sh2ijk);
    gsl_matrix_free(nu_sh2ijk);
    gsl_matrix_free(ijk2cT);
    gsl_matrix_free(sh2cT);
    
#warning call gsl_*_free as needed...
    
    return 0;
}

