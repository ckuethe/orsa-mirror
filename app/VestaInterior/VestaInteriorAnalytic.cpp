#include "VestaInteriorAnalytic.h"

#include <orsa/statistic.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#include "vesta.h"

int main() {
    
    // test specific cases, for debug purposes only!
    // orsa::GlobalRNG::randomSeed = -800402816;
    
    orsa::Debug::instance()->initTimer();
    
    osg::ref_ptr<orsaPDS::RadioScienceGravityFile> pds =
        new orsaPDS::RadioScienceGravityFile("JGDAWN20SIMA.DAT",512,1518);
    
    gsl_vector * pds_coeff    = pds->data->getCoefficientVector();
    gsl_matrix * pds_covm     = pds->data->getCovarianceMatrix();
    gsl_matrix * pds_inv_covm = pds->data->getInverseCovarianceMatrix();
    
    /* if (0) {
       gsl_matrix * identity = gsl_matrix_alloc(pds->data->numberOfCoefficients,pds->data->numberOfCoefficients);
       // gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,pds_covm,pds_inv_covm,0.0,identity);
       gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,pds_inv_covm,pds_covm,0.0,identity);
       for (size_t i=0; i<pds->data->numberOfCoefficients; ++i) {
       for (size_t j=0; j<pds->data->numberOfCoefficients; ++j) {
       ORSA_DEBUG("identity[%03i][%03i] = %+20.9f",i,j,gsl_matrix_get(identity,i,j));
       }
       }
       }
    */
    
    osg::ref_ptr<VestaShape> shape = new VestaShape;
    if (!shape->read("vesta_thomas.dat")) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    const size_t SH_degree = 8; // shperical harmonics degree
    const size_t  T_degree = 6; // chebyshev polinomials degree
    
    const double R0 = pds->data->R0;
#warning which GM value to use? pds->data->GM  OR pds->data->getCoeff("GM") ??
    const double GM = pds->data->GM; 
    
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
                                const size_t z_sh  = (l==0 ? 0 :  pds->data->index(orsaPDS::RadioScienceGravityData::keyC(l,m)));
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

                    
                    
                    // #warning RESTORE
                    const double pq_factor_base = 
                        orsa::power_sign(p+q) *
                        orsa::binomial(l,p).get_d() *
                        orsa::binomial(2*l-2*p,l).get_d() *
                        orsa::binomial(m,2*q).get_d() *
                        orsa::pochhammer(l-m-2*p+1,m);
                    // const double pq_factor_base = 1.0;
                    
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
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, 1.0/int_pow(2.0,l), pq_sh2ijk);
            // gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, 1.0, pq_sh2ijk);
            
            // ORSA_DEBUG("PaulMoment::normalization(%i,%i) = %g",l,m,PaulMoment::normalization(l,m));
            // normalization
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, PaulMoment::normalization(l,m), pq_sh2ijk);


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
                                
                                const size_t z_sh  = pds->data->index(orsaPDS::RadioScienceGravityData::keyS(l,m));
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
                    
                    
                    const double pq_factor_base = 
                        orsa::power_sign(p+q) *
                        orsa::binomial(l,p).get_d() *
                        orsa::binomial(2*l-2*p,l).get_d() *
                        orsa::binomial(m,2*q+1).get_d() *
                        orsa::pochhammer(l-m-2*p+1,m);
	  
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
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, 1.0/int_pow(2.0,l), pq_sh2ijk);
            
            // ORSA_DEBUG("PaulMoment::normalization(%i,%i) = %g",l,m,PaulMoment::normalization(l,m));
            // normalization
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, PaulMoment::normalization(l,m), pq_sh2ijk);
            
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
    
#warning FIX GM entry!!
    
#warning re-include factor 1/R0^l ?
    
#warning re-indlude normalization for Ckm Slm
    
#warning maybe rescale GM entries to 1.0 to have more homogeneous entries?
    
    for (size_t z_sh=0; z_sh<SH_size; ++z_sh) {
        for (size_t z_ijk=0; z_ijk<ijk_size; ++z_ijk) {
            if (gsl_matrix_get(sh2ijk,z_sh,z_ijk)!=0.0) {
                size_t nx,ny,nz;
                CubicChebyshevMassDistribution::triIndex(nx,ny,nz,z_ijk);                
                ORSA_DEBUG("sh2ijk[%03i][%03i] = %+10.6f [%s -> N[%02i][%02i][%02i]]",
                           z_sh,z_ijk,gsl_matrix_get(sh2ijk,z_sh,z_ijk),pds->data->key(z_sh).toStdString().c_str(),nx,ny,nz);
            }
        }
    }
    
    
    // gsl_matrix * cm2bf = gsl_matrix_calloc(ijk_size,ijk_size);
    
    
    CubicChebyshevMassDistribution::CoefficientType coeff;
    CubicChebyshevMassDistribution::resize(coeff,T_degree);
    // dummy coeff, for now
    for (unsigned int i=0; i<=T_degree; ++i) {
        for (unsigned int j=0; j<=T_degree; ++j) {
            for (unsigned int k=0; k<=T_degree; ++k) {
                if (i+j+k<=T_degree) {
                    if ( (i==0) && (j==0) && (k==0) ) {
                        coeff[i][j][k] = 1.0;
                    } else {
                        coeff[i][j][k] = 0.0;
                    }
                }
            }            
        }
    }
    
    osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
        new CubicChebyshevMassDistribution(coeff,R0);
    
#warning use enough points...
    const size_t numSamplePoints = 1000;
    const bool storeSamplePoints = true;
#warning how to manage centerOfMass??
    const double km = orsa::FromUnits(1.0,orsa::Unit::KM);
#warning must SAMPLE on centerOfMass as well! (or is uncertainty too small, and effect negligible?)
    const orsa::Vector centerOfMass = orsa::Vector(0.11*km,-1.15*km,8.50*km);
    
    // at this point the mass distribution is not yet important
    osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape =
        new orsa::RandomPointsInShape(shape,massDistribution,numSamplePoints,storeSamplePoints);
    
    const double volume = orsa::volume(randomPointsInShape);
    
    const double bulkDensity = GM/orsa::Unit::G()/volume;
    const double bulkDensity_gcm3 = orsa::FromUnits(orsa::FromUnits(bulkDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3);
    
    ORSA_DEBUG("bulkDensity coeff: %g",bulkDensity);
    
#warning fix how bulkDensity is used in this code...
    
    gsl_matrix * ijk2cT = gsl_matrix_calloc(ijk_size,T_size);
    
    for (size_t z_cT=0; z_cT<T_size; ++z_cT) {
        
        size_t Tx,Ty,Tz;
        CubicChebyshevMassDistribution::triIndex(Tx,Ty,Tz,z_cT);
        
        for (unsigned int i=0; i<=T_degree; ++i) {
            for (unsigned int j=0; j<=T_degree; ++j) {
                for (unsigned int k=0; k<=T_degree; ++k) {
                    if (i+j+k<=T_degree) {
                        if ( (i==Tx) && (j==Ty) && (k==Tz) ) {
                            coeff[i][j][k] = 1.0;
                        } else {
                            coeff[i][j][k] = 0.0;
                        }
                    }
                }            
            }
        }
        massDistribution = new CubicChebyshevMassDistribution(coeff,R0);
        randomPointsInShape->updateMassDistribution(massDistribution);
        
        for (size_t z_ijk=0; z_ijk<ijk_size; ++z_ijk) {
            
            size_t nx,ny,nz;
            CubicChebyshevMassDistribution::triIndex(nx,ny,nz,z_ijk);
            
            osg::ref_ptr< orsa::Statistic<double> > stat = new orsa::Statistic<double>;
            // osg::ref_ptr< orsa::WeightedStatistic<double> > stat = new orsa::WeightedStatistic<double>;
            
            orsa::Vector v;
            double density;
            
            stat->reset();
            
            randomPointsInShape->reset();
            while (randomPointsInShape->get(v,density)) {
                // const double density = massDistribution->density(v);
                
                // if (density > 0.0) {
                {
                    v -= centerOfMass;
                    // v = shapeToLocal*v;
                    stat->insert(int_pow(v.getX(),nx)*
                                 int_pow(v.getY(),ny)*
                                 int_pow(v.getZ(),nz)*
                                 density);
                    /* stat->insert(int_pow(v.getX(),nx)*
                       int_pow(v.getY(),ny)*
                       int_pow(v.getZ(),nz),
                       density);
                    */
                    // ORSA_DEBUG("density: %g",density);
                }
            }
            
            gsl_matrix_set(ijk2cT,z_ijk,z_cT,stat->average()/int_pow(R0,nx+ny+nz));


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
    
    
    for (size_t z_ijk=0; z_ijk<ijk_size; ++z_ijk) {
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
    
    
    gsl_matrix * sh2cT = gsl_matrix_alloc(SH_size,T_size);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, sh2ijk, ijk2cT, 0.0, sh2cT);
    
    
    for (size_t z_sh=0; z_sh<SH_size; ++z_sh) {
        for (size_t z_cT=0; z_cT<T_size; ++z_cT) {
            // if (gsl_matrix_get(sh2cT,z_sh,z_cT)!=0.0) {
            {             
                size_t Tx,Ty,Tz;
                CubicChebyshevMassDistribution::triIndex(Tx,Ty,Tz,z_cT);
                ORSA_DEBUG("sh2cT[%03i][%03i] = %+9.6f [%s ->  cT[%i][%i][%i]]",
                           z_sh,z_cT,gsl_matrix_get(sh2cT,z_sh,z_cT),pds->data->key(z_sh).toStdString().c_str(),Tx,Ty,Tz);
            }
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
                size_t Tx,Ty,Tz;
                CubicChebyshevMassDistribution::triIndex(Tx,Ty,Tz,s);
                ORSA_DEBUG("uK[%03i][%03i] =%+12.6f (null space base vector for cT[%i][%i][%i])",b,s,gsl_vector_get(uK[b],s),Tx,Ty,Tz);
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
        gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(pds->data->numberOfCoefficients); // workspace for eigenvectors/values
        gsl_vector * eval = gsl_vector_alloc(pds->data->numberOfCoefficients);    // eigenvalues
        gsl_matrix * evec = gsl_matrix_alloc(pds->data->numberOfCoefficients,pds->data->numberOfCoefficients);  // eigenvectors
        //
        gsl_eigen_symmv(pds_covm, eval, evec, w); // NOTE: The diagonal and lower triangular part of A are destroyed during the computation
        //
        double sigma[pds->data->numberOfCoefficients];
        for (size_t i=0; i<pds->data->numberOfCoefficients; ++i) {
            // ORSA_DEBUG("eval[%i] = %g",i,gsl_vector_get(eval,i));
            if (gsl_vector_get(eval,i) == 0.0) {
                ORSA_ERROR("problems with the covariance matrix: null eigenvalue found.");
            }
            sigma[i] = sqrt(fabs(gsl_vector_get(eval,i)));
            // ORSA_DEBUG("sigma[%i] = %g",i,sigma[i]);
        }
        gsl_vector * sampleCoeff_x  = gsl_vector_alloc(pds->data->numberOfCoefficients);
        gsl_vector * sampleCoeff_y  = gsl_vector_alloc(pds->data->numberOfCoefficients); 
        
        for (size_t gen=0; gen<1000; ++gen) {
            
            for (size_t i=0; i<pds->data->numberOfCoefficients; ++i) {
                gsl_vector_set(sampleCoeff_x,i,orsa::GlobalRNG::instance()->rng()->gsl_ran_gaussian(sigma[i]));
            }
            //
            gsl_blas_dgemv(CblasNoTrans,1.0,evec,sampleCoeff_x,0.0,sampleCoeff_y);
            //
            for (size_t i=0; i<pds->data->numberOfCoefficients; ++i) {
                gsl_vector_set(sampleCoeff_y,i,gsl_vector_get(sampleCoeff_y,i)+gsl_vector_get(pds_coeff,i));
            }
            
            for (size_t z_sh=0; z_sh<SH_size; ++z_sh) {
                // gsl_vector_set(sh,z_sh,pds->data->getCoeff(pds->data->key(z_sh)));
                gsl_vector_set(sh,z_sh,gsl_vector_get(sampleCoeff_y,z_sh));
                ORSA_DEBUG("%7s =  %+12.3g [sampled]",pds->data->key(z_sh).toStdString().c_str(),gsl_vector_get(sh,z_sh));
            }

            // solving here!
            gsl_blas_dgemv(CblasNoTrans,1.0,pseudoInvA,sh,0.0,cT);
            
            gsl_vector * cT0 = gsl_vector_calloc(N);
            gsl_vector_memcpy(cT0,cT);

            if (1) {
                // Adaptive Monte Carlo
                osg::ref_ptr<AuxiliaryData> aux = new AuxiliaryData;
                aux->R0 = R0;
                aux->randomPointsInShape = randomPointsInShape;
                aux->bulkDensity = bulkDensity;
                aux->bulkDensity_gcm3 = bulkDensity_gcm3;
                aux->chebyshevDegree = T_degree;
                aux->T_size = T_size;
                aux->cT0 = cT0;
                aux->uK = &uK[0];
                aux->intervalVectorSize = N-M;

                // test
                /* for (size_t z=0; z<N-M; ++z) {
                   for (size_t l=0; l<N; ++l) {
                   ORSA_DEBUG("uK: %g aux->uK: %g",
                   gsl_vector_get(uK[z],l),
                   gsl_vector_get(aux->uK[z],l));
                   }
                   }
                */
                
#warning check this value for DOF!
                size_t chisq_DOF = T_size-SH_size;
                const double chisq_50  = gsl_cdf_chisq_Pinv(0.50,chisq_DOF);
                const double chisq_90  = gsl_cdf_chisq_Pinv(0.90,chisq_DOF);
                const double chisq_95  = gsl_cdf_chisq_Pinv(0.95,chisq_DOF);
                const double chisq_99  = gsl_cdf_chisq_Pinv(0.99,chisq_DOF);
                //
                ORSA_DEBUG("chisq 50\%: %.2f  90\%: %.2f  95\%: %.2f  99\%: %.2f  [DOF=%i]",
                           chisq_50,chisq_90,chisq_95,chisq_99,chisq_DOF);
                
                const double initialThresholdLevel = 1e3; // 100*chisq_99;
                const double targetThresholdLevel  = 0.0; // chisq_99;
                const double maxFactor             = 3.0;
                const double minAdaptiveRange      = -maxFactor;
                const double maxAdaptiveRange      =  maxFactor;
                const double intervalResidualProbability = 1.0e-3; // 1.0e-10; // "1-confidence level" for this interval, different from the chisq-level
                const size_t targetSamples         = 100;
                const size_t maxIter               = 100000;
                
                AdaptiveIntervalVector intervalVector;
                intervalVector.resize(aux->intervalVectorSize);
                for (size_t k=0; k<aux->intervalVectorSize; ++k) {
                    intervalVector[k] =
                        new AdaptiveIntervalType(minAdaptiveRange,
                                                 maxAdaptiveRange,
                                                 intervalResidualProbability,
                                                 initialThresholdLevel,
                                                 targetThresholdLevel,
                                                 targetSamples,
                                                 aux.get());
                }
                
                osg::ref_ptr<AdaptiveMonteCarloType> mc = new AdaptiveMonteCarloType(aux.get());
                
                mc->run(intervalVector,maxIter);
            }
            
            for (unsigned int zuk=0; zuk<1000; ++zuk) {

                if (zuk!=0) {
                    gsl_vector_memcpy(cT,cT0);
                    for (unsigned int b=0; b<(N-M); ++b) {
                        // factor ~ +/- 1000
                        const double maxFactor = 2.0;
                        const double factor = -maxFactor + 2*maxFactor*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
                        for (unsigned int fe=0; fe<N; ++fe) {
                            gsl_vector_set(cT,fe,gsl_vector_get(cT,fe)+factor*gsl_vector_get(uK[b],fe));
                        }
                    }
                }
                
                
                // check for negative density
                {
                    
                    for (unsigned int i=0; i<=T_degree; ++i) {
                        for (unsigned int j=0; j<=T_degree; ++j) {
                            for (unsigned int k=0; k<=T_degree; ++k) {
                                if (i+j+k<=T_degree) {
                                    const size_t index = CubicChebyshevMassDistribution::index(i,j,k);
                                    coeff[i][j][k] = gsl_vector_get(cT,index);
                                }
                            }            
                        }
                    }
                    
                    massDistribution = new CubicChebyshevMassDistribution(coeff,R0);
                    randomPointsInShape->updateMassDistribution(massDistribution);
                    randomPointsInShape->reset();
                    // bool negativeDensity=false;
                    orsa::Vector v;
                    double density;
                    osg::ref_ptr< orsa::Statistic<double> > stat = new orsa::Statistic<double>;
                    size_t numNegativeDensity=0;
                    while (randomPointsInShape->get(v,density)) {
                        stat->insert(density);
                        if (density < 0.0) {
                            // ORSA_DEBUG("negative density...");
                            // negativeDensity=true;
                            ++numNegativeDensity;
                            // don't break, to count all of them
                            // break;
                        }
                    }
                    if (numNegativeDensity!=0) {
                        /* ORSA_DEBUG("negative density [%6u/%u points = %7.3f\%] min: %+6.2f max: %+6.2f avg: %+6.2f [g/cm^3]",
                           numNegativeDensity,randomPointsInShape->size,
                           100.0*(double)numNegativeDensity/(double)randomPointsInShape->size,
                           stat->min(),stat->max(),stat->average());
                        */
                        ORSA_DEBUG("negative density: min: %+6.2f max: %+6.2f avg: %+6.2f [g/cm^3]",
                                   bulkDensity_gcm3*stat->min(),
                                   bulkDensity_gcm3*stat->max(),
                                   bulkDensity_gcm3*stat->average());
                        continue;
                    } else {
                        ORSA_DEBUG("good sample: min: %+6.2f max: %+6.2f avg: %+6.2f [g/cm^3]",
                                   bulkDensity_gcm3*stat->min(),
                                   bulkDensity_gcm3*stat->max(),
                                   bulkDensity_gcm3*stat->average());
                    }
                }
                
                for (size_t z_cT=0; z_cT<T_size; ++z_cT) {
                    size_t Tx,Ty,Tz;
                    CubicChebyshevMassDistribution::triIndex(Tx,Ty,Tz,z_cT);
                    ORSA_DEBUG("cT[%i][%i][%i] = %+12.3f",
                               Tx,Ty,Tz,gsl_vector_get(cT,z_cT));
                }
                
                for (size_t z_cT=0; z_cT<T_size; ++z_cT) {
                    size_t Tx,Ty,Tz;
                    CubicChebyshevMassDistribution::triIndex(Tx,Ty,Tz,z_cT);
                    ORSA_DEBUG("density_coeff_T[%i][%i][%i] = %+12.3f [g/cm^3]",
                               Tx,Ty,Tz,bulkDensity_gcm3*gsl_vector_get(cT,z_cT));
                }
                
                {
                    // one more check, can comment out in production
                    // are the SH coefficient generatead actually those we started with? [sampleCoeff_y]
                    gsl_vector * shReconstructed = gsl_vector_calloc(M); 
                    gsl_blas_dgemv(CblasNoTrans,1.0,sh2cT,cT,0.0,shReconstructed);
                    for (size_t z_sh=0; z_sh<SH_size; ++z_sh) {
                        ORSA_DEBUG("orig: %+12.3g  reconstructed: %+12.3g [%s]",
                                   gsl_vector_get(sh,z_sh),
                                   gsl_vector_get(shReconstructed,z_sh),
                                   pds->data->key(z_sh).toStdString().c_str());
                    }
                }
                
                // output
                for (size_t z_cT=0; z_cT<T_size; ++z_cT) {
                    size_t Tx,Ty,Tz;
                    CubicChebyshevMassDistribution::triIndex(Tx,Ty,Tz,z_cT);
                    gmp_fprintf(stdout,"%u %u %u %+9.6f ",
                                Tx,Ty,Tz,bulkDensity_gcm3*gsl_vector_get(cT,z_cT));
                }
                gmp_fprintf(stdout,"\n");
                fflush(stdout);

                {
                    // raw output of density on file
                    static size_t ID = 0;
                    char filename[1024];
                    sprintf(filename,"density_%06d.dat",ID);
                    ORSA_DEBUG("dumped density to file [%s]",filename);
                    FILE * fp = fopen(filename,"w");
                    if (fp) {
                        orsa::Vector v;
                        double density;
                        randomPointsInShape->reset();
                        while (randomPointsInShape->get(v,density)) {
                            gmp_fprintf(fp,"%+8.3f %+8.3f %+8.3f %+7.3f\n",
                                        orsa::FromUnits(v.getX(),orsa::Unit::KM,-1),
                                        orsa::FromUnits(v.getY(),orsa::Unit::KM,-1),
                                        orsa::FromUnits(v.getZ(),orsa::Unit::KM,-1),
                                        bulkDensity_gcm3*density);
                        }
                        fclose(fp);
                    }
                    ++ID;
                }
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

