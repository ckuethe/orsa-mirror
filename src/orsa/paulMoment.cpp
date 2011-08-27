#include <orsa/paulMoment.h>

#include <orsa/box.h>
#include <orsa/double.h>
// #include <orsa/legendre.h>
#include <orsa/multimin.h>
#include <orsa/statistic.h>
#include <orsa/unit.h>
#include <orsa/util.h>

#include <vector>

#include <gsl/gsl_rng.h>

#include <iostream>

using namespace orsa;

PaulMoment::PaulMoment(const unsigned int n) : osg::Referenced(true), order(n) {
    
    const unsigned int order_plus_one = order+1;
    
    {
        _M.resize(order_plus_one);
        _M_uncertainty.resize(order_plus_one);
        for (unsigned int i=0; i<order_plus_one; ++i) {
            _M[i].resize(order_plus_one-i);
            _M_uncertainty[i].resize(order_plus_one-i);
            for (unsigned int j=0; j<order_plus_one-i; ++j) {
                _M[i][j].resize(order_plus_one-i-j);
                _M_uncertainty[i][j].resize(order_plus_one-i-j);
            }
        }
    }
}

double PaulMoment::M (const int i,
                      const int j, 
                      const int k) const {
    if ( (i<0) || (j<0) || (k<0) || (i+j+k > (int)order) ) { 
        // ORSA_ERROR("index out of range [i:%i;j:%i;k:%i]",i,j,k);
        return 0; 
    } else {
        return _M[i][j][k];
    }
}

double PaulMoment::M_uncertainty (const int i,
                                  const int j, 
                                  const int k) const {
    if ( (i<0) || (j<0) || (k<0) || (i+j+k > (int)order) ) { 
        // ORSA_ERROR("index out of range [i:%i;j:%i;k:%i]",i,j,k);
        return 0;
    } else {
        return _M_uncertainty[i][j][k];
    }
}

void PaulMoment::setM (const double & val,
                       const int i, 
                       const int j, 
                       const int k) {
    if ( (i<0) || (j<0) || (k<0) || (i+j+k > (int)order) ) {
        ORSA_ERROR("index out of range [i:%i;j:%i;k:%i]",i,j,k);
        return;
    } else {
        _M[i][j][k] = val;
    }
}

void PaulMoment::setM_uncertainty (const double & val,
                                   const int i, 
                                   const int j, 
                                   const int k) {
    if ( (i<0) || (j<0) || (k<0) || (i+j+k > (int)order) ) { 
        ORSA_ERROR("index out of range [i:%i;j:%i;k:%i]",i,j,k);
        return;
    } else {
        _M_uncertainty[i][j][k] = val;
    }
}

/***/

const orsa::triIndex_mpq orsa::conversionCoefficients_C_integral(const size_t & l_ask, const size_t & m_ask) {
    static std::deque< std::deque<orsa::triIndex_mpq> > coeff;
    if (coeff.size() > l_ask) {
        if (coeff[l_ask].size() > m_ask) {
            return coeff[l_ask][m_ask];
        }        
    }
    const size_t old_l_size = coeff.size();
    coeff.resize(l_ask+1);
    for (int l=old_l_size; l<=(int)l_ask; ++l) {
        coeff[l].resize(l+1);
        const mpz_class pow_2_l = orsa::int_pow(mpz_class(2),l);
        for (int m=0; m<=l; ++m) {
            
            coeff[l][m].resize(l+1);
            for (int ti=0; ti<=l; ++ti) {
                coeff[l][m][ti].resize(l+1-ti);
                for (int tj=0; tj<=l-ti; ++tj) {
                    coeff[l][m][ti][tj].resize(l+1-ti-tj);
                }
            }
            
            triIndex_mpq pq_factor;
            pq_factor.resize(l+1);
            for (int ti=0; ti<=l; ++ti) {
                pq_factor[ti].resize(l+1-ti);
                for (int tj=0; tj<=l-ti; ++tj) {
                    pq_factor[ti][tj].resize(l+1-ti-tj);
                }
            }
            
            for (int ti=0; ti<=l; ++ti) {
                for (int tj=0; tj<=l-ti; ++tj) {
                    for (int tk=0; tk<=l-ti-tj; ++tk) {
                        pq_factor[ti][tj][tk] = 0;
                    }   
                }
            }
            
            // integer division in limits
            for (int p=0;p<=(l/2);++p) {
                for (int q=0;q<=(m/2);++q) {

                    triIndex_mpq nu_factor;
                    nu_factor.resize(l+1);
                    for (int ti=0; ti<=l; ++ti) {
                        nu_factor[ti].resize(l+1-ti);
                        for (int tj=0; tj<=l-ti; ++tj) {
                            nu_factor[ti][tj].resize(l+1-ti-tj);
                        }
                    }
                    
                    for (int ti=0; ti<=l; ++ti) {
                        for (int tj=0; tj<=l-ti; ++tj) {
                            for (int tk=0; tk<=l-ti-tj; ++tk) {
                                nu_factor[ti][tj][tk] = 0;
                            }   
                        }
                    }
                    
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
                                
                                const mpq_class nu_factor_base(orsa::factorial(p),orsa::factorial(nu_x)*orsa::factorial(nu_y)*orsa::factorial(p-nu_x-nu_y));
                                
                                /* ORSA_DEBUG("nu_factor_base[%i][%i][%i] += %Zi/%Zi = %i!/(%i!%i!%i!)",
                                   M_i, M_j, M_k,
                                   nu_factor_base.get_num().get_mpz_t(),
                                   nu_factor_base.get_den().get_mpz_t(),
                                   p,nu_x,nu_y,p-nu_x-nu_y);
                                */
                                
                                nu_factor[M_i][M_j][M_k] += nu_factor_base;
                            }
                        }
                    }
                    
                    const mpz_class pq_factor_base = 
                        orsa::power_sign(p+q) *
                        orsa::binomial(l,p) *
                        orsa::binomial(2*l-2*p,l) *
                        orsa::binomial(m,2*q) *
                        orsa::pochhammer(mpz_class(l-m-2*p+1),m);
                    
                    for (int ti=0; ti<=l; ++ti) {
                        for (int tj=0; tj<=l-ti; ++tj) {
                            for (int tk=0; tk<=l-ti-tj; ++tk) {
                                nu_factor[ti][tj][tk] *= pq_factor_base;
                            }   
                        }
                    }
                    
                    for (int ti=0; ti<=l; ++ti) {
                        for (int tj=0; tj<=l-ti; ++tj) {
                            for (int tk=0; tk<=l-ti-tj; ++tk) {
                                pq_factor[ti][tj][tk] += nu_factor[ti][tj][tk];
                            }   
                        }
                    }
                }
            }
            
            for (int ti=0; ti<=l; ++ti) {
                for (int tj=0; tj<=l-ti; ++tj) {
                    for (int tk=0; tk<=l-ti-tj; ++tk) {
                        // divide by 2^l
                        pq_factor[ti][tj][tk] /= pow_2_l;
                    }   
                }
            }
            
            for (int ti=0; ti<=l; ++ti) {
                for (int tj=0; tj<=l-ti; ++tj) {
                    for (int tk=0; tk<=l-ti-tj; ++tk) {
                        pq_factor[ti][tj][tk].canonicalize();
                        coeff[l][m][ti][tj][tk] = pq_factor[ti][tj][tk];
                    }
                }
            }
        }
    }
    return coeff[l_ask][m_ask];
}

/***/

void orsa::convert(std::vector< std::vector<mpf_class> > & C,
                   std::vector< std::vector<mpf_class> > & S,
                   std::vector< std::vector<mpf_class> > & norm_C,
                   std::vector< std::vector<mpf_class> > & norm_S,
                   std::vector<mpf_class> & J,
                   const PaulMoment * const pm,
                   const double     & R0,
                   const bool         verbose) {
  
    if (verbose) {
        ORSA_DEBUG("R0: %f [km]",orsa::FromUnits(R0,orsa::Unit::KM,-1));
    }
  
    const unsigned int order = pm->order;
  
    // resize vectors
    C.resize(order+1);
    S.resize(order+1);
    norm_C.resize(order+1);
    norm_S.resize(order+1);
    //
    for (unsigned int l=0; l<=order; ++l) {
        C[l].resize(l+1);
        S[l].resize(l+1);
        norm_C[l].resize(l+1);
        norm_S[l].resize(l+1);
    }
    //
    J.resize(order+1);
  
    // C_lm coefficients
    //
    for (int l=0; l<=(int)order; ++l) {
        for (int m=0; m<=l; ++m) {
      
            mpf_class pq_factor=0;
            mpf_class pq_factor_uncertainty=0;
            //
            // integer division in limits
            for (int p=0;p<=(l/2);++p) {
                for (int q=0;q<=(m/2);++q) {
	  
                    mpf_class nu_factor=0;
                    mpf_class nu_factor_uncertainty=0;
                    //
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
		
                                const mpf_class nu_factor_base =
                                    (orsa::factorial(p).get_d() /(orsa::factorial(nu_x).get_d()*orsa::factorial(nu_y).get_d()*orsa::factorial(p-nu_x-nu_y).get_d()));
		
                                nu_factor += nu_factor_base * pm->M(M_i,M_j,M_k);
                                nu_factor_uncertainty += nu_factor_base * pm->M_uncertainty(M_i,M_j,M_k);
                            }
                        }
                    }
                    // need this because using M(i,j,k) instead of N(i,j,k)
                    nu_factor /= int_pow(R0,l);
                    nu_factor_uncertainty /= int_pow(R0,l);
	  
                    const mpf_class pq_factor_base = 
                        orsa::power_sign(p+q) *
                        orsa::binomial(l,p).get_d() *
                        orsa::binomial(2*l-2*p,l).get_d() *
                        orsa::binomial(m,2*q).get_d() *
                        orsa::pochhammer(mpz_class(l-m-2*p+1),m);
	  
                    pq_factor += 
                        pq_factor_base *
                        nu_factor;
	  
                    pq_factor_uncertainty += 
                        pq_factor_base *
                        nu_factor_uncertainty;
                }
            }
            
            pq_factor /= int_pow(2.0,l);
            pq_factor_uncertainty /= int_pow(2.0,l);
            
            const mpf_class C_lm = orsa::normalization_integralToSphericalHarmonics(l,m)*pq_factor;
            const mpf_class C_lm_uncertainty = orsa::normalization_integralToSphericalHarmonics(l,m)*abs(pq_factor_uncertainty);
            //
            const mpf_class norm_C_lm = C_lm*orsa::normalization_sphericalHarmonicsToNormalizedSphericalHarmonics(l,m);
            const mpf_class norm_C_lm_uncertainty = abs(C_lm_uncertainty*orsa::normalization_sphericalHarmonicsToNormalizedSphericalHarmonics(l,m));
            
            if (verbose) {
                ORSA_DEBUG("     C[%i][%i] = %+16.12Ff +/- %16.12Ff",
                           l,m,     C_lm.get_mpf_t(),     C_lm_uncertainty.get_mpf_t());
                ORSA_DEBUG("norm_C[%i][%i] = %+16.12Ff +/- %16.12Ff   norm: %Ff",
                           l,m,norm_C_lm.get_mpf_t(),norm_C_lm_uncertainty.get_mpf_t(),
                           orsa::normalization_sphericalHarmonicsToNormalizedSphericalHarmonics(l,m).get_mpf_t());
            }
            
            C[l][m]      = C_lm;
            norm_C[l][m] = norm_C_lm;
      
            // J_l is minus C_l0, where C_l0 is not normalized
            if (l>=2) {
                if (m==0) {
                    const mpf_class J_l = -C_lm;
                    const mpf_class J_l_uncertainty = -C_lm_uncertainty;
                    if (verbose) {
                        ORSA_DEBUG("J_%i = %+16.12Ff +/- %16.12Ff",l,J_l.get_mpf_t(),J_l_uncertainty.get_mpf_t());
                    }
                    //
                    J[l] = J_l;
                }	
            } else {
                if (m==0) {
                    J[l] = 0.0;
                }
            }
      
        }
    }
  
    // S_lm coefficients
    //
    for (int l=0; l<=(int)order; ++l) {
        for (int m=1; m<=l; ++m) {
      
            mpf_class pq_factor=0;
            mpf_class pq_factor_uncertainty=0;
            //
            // integer division in limits
            for (int p=0;p<=(l/2);++p) {
                for (int q=0;q<=((m-1)/2);++q) {
	  
                    mpf_class nu_factor=0;
                    mpf_class nu_factor_uncertainty=0;
                    //
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
		
                                const mpf_class nu_factor_base =
                                    (orsa::factorial(p).get_d() /(orsa::factorial(nu_x).get_d()*orsa::factorial(nu_y).get_d()*orsa::factorial(p-nu_x-nu_y).get_d()));
		
                                nu_factor += nu_factor_base * pm->M(M_i,M_j,M_k);
                                nu_factor_uncertainty += nu_factor_base * pm->M_uncertainty(M_i,M_j,M_k);
                            }
                        }
                    }
                    // need this because using M(i,j,k) instead of N(i,j,k)
                    nu_factor /= int_pow(R0,l);
                    nu_factor_uncertainty /= int_pow(R0,l);
	  
                    const mpf_class pq_factor_base = 
                        orsa::power_sign(p+q) *
                        orsa::binomial(l,p).get_d() *
                        orsa::binomial(2*l-2*p,l).get_d() *
                        orsa::binomial(m,2*q+1).get_d() *
                        orsa::pochhammer(mpz_class(l-m-2*p+1),m);
                    
                    pq_factor += 
                        pq_factor_base *
                        nu_factor;
	  
                    pq_factor_uncertainty += 
                        pq_factor_base *
                        nu_factor_uncertainty;
                }
            }
            //
            pq_factor /= int_pow(2.0,l);
            pq_factor_uncertainty /= int_pow(2.0,l);
            //
            const mpf_class S_lm = orsa::normalization_integralToSphericalHarmonics(l,m)*pq_factor;
            const mpf_class S_lm_uncertainty = orsa::normalization_integralToSphericalHarmonics(l,m)*abs(pq_factor_uncertainty);
            //      
            const mpf_class norm_S_lm = S_lm*orsa::normalization_sphericalHarmonicsToNormalizedSphericalHarmonics(l,m);
            const mpf_class norm_S_lm_uncertainty = abs(S_lm_uncertainty*orsa::normalization_sphericalHarmonicsToNormalizedSphericalHarmonics(l,m));
      
            if (verbose) {
                ORSA_DEBUG("     S[%i][%i] = %+16.12Ff +/- %16.12Ff",
                           l,m,     S_lm.get_mpf_t(),     S_lm_uncertainty.get_mpf_t());
                ORSA_DEBUG("norm_S[%i][%i] = %+16.12Ff +/- %16.12Ff   norm: %Ff",
                           l,m,norm_S_lm.get_mpf_t(),norm_S_lm_uncertainty.get_mpf_t(),
                           orsa::normalization_sphericalHarmonicsToNormalizedSphericalHarmonics(l,m).get_mpf_t());
            }
      
            S[l][m]      = S_lm;
            norm_S[l][m] = norm_S_lm;
        }
    }
  
}

/*** Use Multimin code to solve the inverse problem ***/

class PaulMomentsSolveMultimin : public orsa::Multimin {
public:
    PaulMomentsSolveMultimin(const unsigned int focusOrder_in, 
                             const std::vector< std::vector<mpf_class> > & norm_C_in,
                             const std::vector< std::vector<mpf_class> > & norm_S_in,
                             const double & R0_in) : 
        orsa::Multimin(), 
        focusOrder(focusOrder_in),
        norm_C(norm_C_in),
        norm_S(norm_S_in),
        R0(R0_in) { }
protected:
    const unsigned int focusOrder;
    const std::vector< std::vector<mpf_class> > & norm_C;
    const std::vector< std::vector<mpf_class> > & norm_S;
    const double R0;
public:
    double fun(const orsa::MultiminParameters * par) const {
        osg::ref_ptr<PaulMoment> local_pm = new PaulMoment(focusOrder);
        char parName[1024];
        for (unsigned int i=0; i<=focusOrder; ++i) {
            for (unsigned int j=0; j<=focusOrder; ++j) {
                for (unsigned int k=0; k<=focusOrder; ++k) {
                    if (i+j+k==focusOrder) {
                        sprintf(parName,"M_%i_%i_%i",i,j,k);
                        local_pm->setM(par->get(parName),i,j,k);
                    }
                }
            }
        }
    
        // new C,S values
        std::vector< std::vector<mpf_class> > local_C, local_S, local_norm_C, local_norm_S;
        std::vector<mpf_class> local_J;
        convert(local_C, local_S, local_norm_C, local_norm_S, local_J,
                local_pm.get(),
                R0);
    
        mpf_class retVal=0.0;
        {
            const unsigned int l=focusOrder;
            for (unsigned int m=0; m<=l; ++m) {
                retVal += orsa::square(mpf_class(local_norm_C[l][m]-norm_C[l][m]));
                if (m!=0) {
                    retVal += orsa::square(mpf_class(local_norm_S[l][m]-norm_S[l][m]));	  
                }
            }
        }
    
        return retVal.get_d();
    }
};

bool orsa::solve(PaulMoment * pm,
                 const std::vector< std::vector<mpf_class> > & norm_C,
                 const std::vector< std::vector<mpf_class> > & norm_S,
                 const double     & R0) {
  
    const unsigned int order = pm->order;
  
    for (unsigned int focusOrder=0; focusOrder<=order; ++focusOrder) {
    
        osg::ref_ptr<orsa::MultiminParameters> par = new orsa::MultiminParameters;
        //
        {
            char parName[1024];
            // for (unsigned int sum=0; sum<=order; ++sum) {
            for (unsigned int i=0; i<=focusOrder; ++i) {
                for (unsigned int j=0; j<=focusOrder; ++j) {
                    for (unsigned int k=0; k<=focusOrder; ++k) {
                        if (i+j+k==focusOrder) {
                            sprintf(parName,"M_%i_%i_%i",i,j,k);
                            // par->insert(parName,0.01*int_pow(R0,i+j+k),1e16*(int_pow(R0,i+j+k)));
                            par->insert(parName,0.01*int_pow(R0,focusOrder),int_pow(R0,focusOrder));
                        }
                    }
                }
            }
            // }
        }
        
        osg::ref_ptr<PaulMomentsSolveMultimin> multimin = 
            new PaulMomentsSolveMultimin(focusOrder,norm_C,norm_S,R0);
        //
        multimin->setMultiminParameters(par.get());
        //
        multimin->run_nmsimplex(1024*1024,1.0e-15*int_pow(R0,focusOrder));
        // multimin->run_conjugate_fr(1024*1024,1.0,1e-6,1.0e-15*int_pow(R0,focusOrder));
    
        // save results
        {
            char parName[1024];
            for (unsigned int i=0; i<=focusOrder; ++i) {
                for (unsigned int j=0; j<=focusOrder; ++j) {
                    for (unsigned int k=0; k<=focusOrder; ++k) {
                        if (i+j+k==focusOrder) { 
                            sprintf(parName,"M_%i_%i_%i",i,j,k);
                            pm->setM(par->get(parName),i,j,k);
                            pm->setM_uncertainty(0,i,j,k);
                        }
                    }
                }
            }
        }
    
    }
  
    return true;
}

static mpz_class EllipsoidExpansion_product_utility(const unsigned int n) {
    mpz_class product = 1;
    for (unsigned int k=1; k<=n; ++k) {
        product *= (2*k-1);
    }
    return product;
}

void orsa::EllipsoidExpansion(PaulMoment   * pm,
                              const double & a,
                              const double & b,
                              const double & c) {
  
    const unsigned int order = pm->order;
  
    for (unsigned int focusOrder=0; focusOrder<=order; ++focusOrder) {
        for (unsigned int i=0; i<=focusOrder; ++i) {
            for (unsigned int j=0; j<=focusOrder; ++j) {
                for (unsigned int k=0; k<=focusOrder; ++k) {
                    if (i+j+k==focusOrder) {
                        if (i%2==1) continue;
                        if (j%2==1) continue;
                        if (k%2==1) continue;
                        const mpz_class factor_i   = EllipsoidExpansion_product_utility(i/2);
                        const mpz_class factor_j   = EllipsoidExpansion_product_utility(j/2);
                        const mpz_class factor_k   = EllipsoidExpansion_product_utility(k/2);
                        const mpz_class factor_ijk = EllipsoidExpansion_product_utility((i+j+k)/2+2);
                        //
                        mpq_class factor(3*(factor_i*factor_j*factor_k),factor_ijk);
                        factor.canonicalize();
                        //
                        const double M_ijk         = factor.get_d()*orsa::int_pow(a,i)*orsa::int_pow(b,j)*orsa::int_pow(c,k);
                        pm->setM(M_ijk,i,j,k);
                        ORSA_DEBUG("ijk: %i %i %i   factor: %Zi/%Zi  M: %g",
                                   i,j,k,
                                   factor.get_num().get_mpz_t(),factor.get_den().get_mpz_t(),
                                   M_ijk); 
                    }
                }
            }
        }
    }
}
