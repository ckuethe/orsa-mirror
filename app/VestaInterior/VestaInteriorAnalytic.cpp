#include "VestaInteriorAnalytic.h"

#include <orsa/statistic.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>

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
    
    const size_t SH_degree = 4; // shperical harmonics degree
    const size_t  T_degree = 2; // chebyshev polinomials degree

    const double R0 = pds->data->R0;
    
    // sh size: (l+1)^2 +1; +1 due to GM factor; -4 because C00, C10, C11, S11 are missing
    const size_t  SH_size = (SH_degree+1)*(SH_degree+1)+1-4;
    const size_t ijk_size = CubicChebyshevMassDistribution::totalSize(SH_degree);
    const size_t   T_size = CubicChebyshevMassDistribution::totalSize( T_degree);
    
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

                                const size_t z_sh  = pds->data->index(orsaPDS::RadioScienceGravityData::keyC(l,m));
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
            
            ORSA_DEBUG("PaulMoment::normalization(%i,%i) = %g",l,m,PaulMoment::normalization(l,m));
            // normalization
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, identity_sh, pq_sh2ijk, PaulMoment::normalization(l,m), pq_sh2ijk);
            
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
            
            ORSA_DEBUG("PaulMoment::normalization(%i,%i) = %g",l,m,PaulMoment::normalization(l,m));
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
    
#warning check for negative densities!
    
#warning re-include factor 1/R0^l ?
    
#warning re-indlude normalization for Ckm Slm
    
    for (size_t z_sh=0; z_sh<SH_size; ++z_sh) {
        for (size_t z_ijk=0; z_ijk<ijk_size; ++z_ijk) {
            if (gsl_matrix_get(sh2ijk,z_sh,z_ijk)!=0.0) {
                size_t nx,ny,nz;
                CubicChebyshevMassDistribution::triIndex(nx,ny,nz,z_ijk);                
                ORSA_DEBUG("sh2ijk[%02i][%02i] = %+12.3g   [%s -> N[%02i][%02i][%02i]]",
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
    const size_t numSamplePoints = 10000;
    const bool storeSamplePoints = true;
#warning how to manage centerOfMass??
    const double km = orsa::FromUnits(1.0,orsa::Unit::KM);                                      
    const orsa::Vector centerOfMass = orsa::Vector(0.0*km,0.0*km,5.0*km);
    
    // at this point the mass distribution is not yet important
    osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape =
        new orsa::RandomPointsInShape(shape,massDistribution,numSamplePoints,storeSamplePoints);
    
    const double volume = orsa::volume(randomPointsInShape);
    
    gsl_matrix * ijk2cT = gsl_matrix_calloc(ijk_size,T_size);
    
    for (size_t z_ijk=0; z_ijk<ijk_size; ++z_ijk) {
        
        size_t nx,ny,nz;
        CubicChebyshevMassDistribution::triIndex(nx,ny,nz,z_ijk);
        
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

            osg::ref_ptr< orsa::WeightedStatistic<double> > stat = new orsa::WeightedStatistic<double>;
            
            orsa::Vector v;
            double density;
            
            stat->reset();
            
            randomPointsInShape->reset();
            while (randomPointsInShape->get(v,density)) {
                // const double density = massDistribution->density(v);
                if (density > 0) {
                    v -= centerOfMass;
                    // v = shapeToLocal*v;
                    stat->insert(int_pow(v.getX(),nx)*
                                 int_pow(v.getY(),ny)*
                                 int_pow(v.getZ(),nz),
                                 density);
                }
            }
            
            gsl_matrix_set(ijk2cT,z_ijk,z_cT,stat->average());
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
                ORSA_DEBUG("ijk2cT[%02i][%02i] = %+12.3g   [N[%02i][%02i][%02i] -> c*Tx[%i][%i][%i]]",
                           z_ijk,z_cT,gsl_matrix_get(ijk2cT,z_ijk,z_cT),nx,ny,nz,Tx,Ty,Tz);
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
    
#warning call gsl_*_free as needed...
    
    return 0;
}

