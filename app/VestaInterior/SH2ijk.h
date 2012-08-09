#ifndef _SH2IJK_H_
#define _SH2IJK_H_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/double.h>
#include <orsa/shape.h>
#include <orsa/unit.h>
#include <orsa/util.h>

#include <qd/dd_real.h>
#include <qd/qd_real.h>

// SQLite3
#include "sqlite3.h"

// #include <gsl/gsl_sf_gamma.h>

#warning should change the _simplex_ part of the name to _ijk_ or something like that
inline std::string getSqliteDBFileName_SH(const std::string & inputFile,
                                          const double & R0) {
    char line[1024];
    sprintf(line,"%s_simplex_%gkm.sqlite",inputFile.c_str(),orsa::FromUnits(R0,orsa::Unit::KM,-1));
    return line;
}

// from src/orsa/double.cpp
template <class T> T orsa::int_pow(const T & x,
                                   const int & p) {
    if (p ==  2) return x*x;
    if (p ==  1) return x;
    if (p ==  0) return 1;
    if (p == -1) return 1/x;
    T _pow = x;
    const int max_k = abs(p);
    for (int k=1; k<max_k; ++k) {
        _pow *= x;
    }
    if (p < 0) _pow = T(1)/_pow;
    return _pow;
}

template dd_real orsa::int_pow(const dd_real & z,
                               const int & p);

typedef std::vector< std::vector<double> > SHcoeff;

inline void writeSH(const SHcoeff & norm_A,
                    const SHcoeff & norm_B,
                    const std::string & outputFile,
                    const orsa::Unit::LengthUnit & outputFileLengthScale = orsa::Unit::METER) {
    ORSA_DEBUG("need code here!");
}

inline void readSH(SHcoeff & norm_A,
                   SHcoeff & norm_B,
                   const std::string & inputFile,
                   const orsa::Unit::LengthUnit & inputFileLengthScale = orsa::Unit::METER) {
    FILE * fp = fopen(inputFile.c_str(),"r");
    if (!fp) {
        ORSA_DEBUG("cannot open file [%s]",inputFile.c_str());
        return;
    }
    
    norm_A.clear();
    norm_B.clear();

    char line[4096];
    char s_l[4096], s_m[4096], s_C[4096], s_S[4096];
    while (fgets(line,4096,fp)) {
        int num_good = sscanf(line,
                              "%s %s %s %s",
                              s_l, s_m, s_C, s_S);
        // ORSA_DEBUG("[%s] [%s] [%s] [%s] (num_good=%i)",s_l, s_m, s_C, s_S, num_good);
        if ( (num_good == 3) ||
             (num_good == 4) ) {
            // test for good lines
            if (!isdigit(s_l[0])) continue;

            const int l=atoi(s_l);
            const int m=atoi(s_m);
            
            if (norm_A.size()<=l) {
                norm_A.resize(l+1);
                norm_B.resize(l+1);
                for (size_t j=0; j<=l; ++j) {
                    norm_A[j].resize(j+1);
                    norm_B[j].resize(j+1);
                }
            }
            
            norm_A[l][m] = orsa::FromUnits(atof(s_C),inputFileLengthScale);
            if ( (num_good==4) && (m>0) ) {
                norm_B[l][m] = orsa::FromUnits(atof(s_S),inputFileLengthScale);
            }   
        }
    }
    
    /* for (size_t l=0; l<norm_A.size(); ++l) {
       for (size_t m=0; m<=l; ++m) {
       ORSA_DEBUG("norm_A[%i][%i] = %g",l,m,norm_A[l][m]);
       ORSA_DEBUG("norm_B[%i][%i] = %g",l,m,norm_B[l][m]);
       }
       }
    */
    
    fclose(fp);
}

template <typename T> class SHIntegration : public osg::Referenced {
public:
    SHIntegration(const SHcoeff & nA,
                  const SHcoeff & nB,
                  const T & R0_,
                  const std::string & SQLiteDBFileName) :
        osg::Referenced(),
        // N(3),
        norm_A(nA),
        norm_B(nB),
        // triShape(s),
        // R0(R0_),
        oneOverR0(1.0/R0_) {
        
        int rc = sqlite3_open(SQLiteDBFileName.c_str(),&db);
        //
        if (rc) {
            fprintf(stderr,"Can't open db: %s\n",sqlite3_errmsg(db));
            sqlite3_close(db);
        }
        
        {
            // get list of tables, to see if some results have been obtained already
            // if no table is present, create it
            char **result;
            int nrows, ncols;
            char * zErr;
            std::string sql = "SELECT name FROM sqlite_master";
            rc = sqlite3_get_table(db,sql.c_str(),&result,&nrows,&ncols,&zErr);
            //
            if (rc != SQLITE_OK) {
                if (zErr != NULL) {
                    fprintf(stderr,"SQL error: %s\n",zErr);
                    sqlite3_free(zErr);
                }
            }
            // ORSA_DEBUG("nrows: %i  ncols: %i",nrows, ncols);
            //
            /* for (int i=0; i<nrows; ++i) {
               for (int j=0; j<ncols; ++j) {
               // i=0 is the header
               const int index = (i+1)*ncols+j;
               ORSA_DEBUG("result[%i] = %s",index, result[index]);
               }
               }
            */
            //
            bool createTable=true;
            //
            // if nrows == number of tables, don't need to create them
            if (nrows==1) {
                createTable=false;
            }
            //
            sqlite3_free_table(result);
            
            if (createTable) {
                // create results table
                sql = "CREATE TABLE simplex(id INTEGER PRIMARY KEY, nx INTEGER, ny INTEGER, nz INTEGER, integral REAL)";
                rc = sqlite3_exec(db,sql.c_str(),NULL,NULL,&zErr);
                //
                if (rc != SQLITE_OK) {
                    if (zErr != NULL) {
                        fprintf(stderr,"SQL error: %s\n",zErr);
                        sqlite3_free(zErr);
                    }
                }
            }
        }
    }
protected:
    virtual ~SHIntegration() {
        sqlite3_close(db);
    }

public:
    void reserve(const size_t maxDegree) {
        const size_t maxIndex = std::max(getIndex(maxDegree,0,0),
                                         std::max(getIndex(0,maxDegree,0),
                                                  getIndex(0,0,maxDegree)));
        val.reserve(maxIndex+1);
        // val_vol_sum_fun.reserve(maxIndex+1); 
        // pochhammer_Np1.reserve(maxDegree+1);
    }
protected:
    /* class SHInternals {
       public:
       std::vector<orsa::Vector> simplexVertexVector;
       orsa::Cache<T> volume;
       };
    */
protected:
    // const size_t N; // N = 3 = dimension of space = number of vertexes in n-dim simplex + 1
    // osg::ref_ptr<const orsa::TriShape> triShape;
    const SHcoeff & norm_A;
    const SHcoeff & norm_B;
    // const double R0;
    const T oneOverR0;
    mutable std::vector< orsa::Cache<double> > val; // integral value
    // mutable std::vector< SHInternals > aux;
    // val_vol_sum_fun is the sum over all simplexes of the funciton of given q times the volume of each simplex
    // mutable std::vector< std::vector< orsa::Cache<T> > > val_vol_sum_fun;
    // mutable std::vector< orsa::Cache<mpz_class> > pochhammer_Np1; // pochhammer_Np1[deg] = pochhammer(N+1,deg)
protected:
    // needed to work with SQLite database
    sqlite3 * db;
protected:
    // utility functions to allow the use of templates
    static dd_real mpzToDD(const mpz_class & z) {
        // return dd_real(z.get_d());
        char * str = (char *)malloc(mpz_sizeinbase(z.get_mpz_t(),10)+2);
        mpz_get_str(str,10,z.get_mpz_t());
        dd_real x(str);
        // ORSA_DEBUG("z: [%Zi]  STR: [%s]   dd: %g",z.get_mpz_t(),str,::to_double(x));
        free(str);
        return x;
    }
    static qd_real mpzToQD(const mpz_class & z) {
        char * str = (char *)malloc(mpz_sizeinbase(z.get_mpz_t(),10)+2);
        mpz_get_str(str,10,z.get_mpz_t());
        qd_real x(str);
        // ORSA_DEBUG("z: [%Zi]  STR: [%s]   qd: %g",z.get_mpz_t(),str,::to_double(x));
        free(str);
        return x;
    }
    
    template <class U> T aux_01(const int & sign, const mpz_class & binomial, const U & val) const;
    double aux_01(const int & sign, const mpz_class & binomial, const double & val) const { return sign*binomial.get_d()*val; }   
    mpf_class aux_01(const int & sign, const mpz_class & binomial, const mpf_class & val) const { return sign*binomial*val; }   
    dd_real aux_01(const int & sign, const mpz_class & binomial, const dd_real & val) const { return sign*mpzToDD(binomial)*val; }  
    qd_real aux_01(const int & sign, const mpz_class & binomial, const qd_real & val) const { return sign*mpzToQD(binomial)*val; }   
    //
    template <class U> T aux_02(const U & x, const mpz_class & pochhammer) const;
    double aux_02 (const double & x, const mpz_class & pochhammer) const {
        return x/pochhammer.get_d();
    }
    double aux_02 (const mpf_class & x, const mpz_class & pochhammer) const {
        return mpf_class(x/pochhammer).get_d();
    }
    double aux_02 (const dd_real & x, const mpz_class & pochhammer) const {
        dd_real y = x / mpzToDD(pochhammer);
        return ::to_double(y);
    }
    double aux_02 (const qd_real & x, const mpz_class & pochhammer) const {
        qd_real y = x / mpzToQD(pochhammer);
        return ::to_double(y);
    }
    
    // return Gamma(n) = (n-1)!
    dd_real aux_gamma_n(const mpz_class & n) const {
        // ORSA_DEBUG("%g",::to_double(mpzToDD(orsa::factorial(n-1))));
        return mpzToDD(orsa::factorial(n-1));
    }
    
    // return Gamma(n+1/2)
    dd_real aux_gamma_n_plus_half(const mpz_class & n) const {
        dd_real result=1;
        result *= mpzToDD(orsa::factorial(2*n));
        result /= mpzToDD(orsa::int_pow((mpz_class)4,n)*orsa::factorial(n));
        result *= sqrt(dd_real::_pi);
        return result;
    }
    
    // Gamma(n/2)
    dd_real aux_gamma_half_n(const mpz_class & n) const {
        // ORSA_DEBUG("n: %Zi",n.get_mpz_t());
        if (n%2==0) return aux_gamma_n(n/2);
        else return aux_gamma_n_plus_half(n/2);
    }
    
    // Gamma((c+1)/2)*Gamma((s+1)/2)/Gamma((c+s+2)/2)
    dd_real triple_factorial(const mpz_class & c, const mpz_class & s) const {
        // ORSA_DEBUG("c: %Zi  s: %Zi",c.get_mpz_t(),s.get_mpz_t());
        return (aux_gamma_half_n(c+1)*aux_gamma_half_n(s+1)/aux_gamma_half_n(c+s+2));
    }
    
    // integral between 0 and k*pi of cos^c(x) sin^s(x) dx
    dd_real integral_csk(const int & c, const int & s, const int & k) const {
        // ORSA_DEBUG("csk: %i %i %i",c,s,k);
        if ((c<0) || (s<0) || (k<0)) return 0; // tested only for natural numbers values
        const bool c_even = (c%2==0);
        if (c_even) {
            const bool s_even = (s%2==0);
            if (s_even) {
                return k*triple_factorial(c,s);
            } else {
                const bool k_even = (k%2==0);
                if (k_even) {
                    return 0.0;
                } else {
                    return triple_factorial(c,s);
                }
            }
        } else {
            return 0.0;
        }
    }
    
    // norm_coeff = normalization_factor * coeff
    static dd_real normalization_factor(const int & l,
                                        const int & m) {
        /* static std::vector< std::vector< orsa::Cache<dd_real> > > stored;
           if (stored.size() > l) {
           // if (stored[l].size() > m) {
           if (stored[l][m].isSet()) {
           return stored[l][m];
           }
           // }
           } else {
           stored.resize(l+1);
           for (int j=0; j<=l; ++j) {
           stored[j].resize(j+1);
           }
           }
        */
        
        // return mpfToDD(orsa::normalization_sphericalHarmonicsToNormalizedSphericalHarmonics(l,m));
        // write it explicitely to avoid losing digits...
        const dd_real val =
            sqrt(mpzToDD(orsa::factorial(l+m)) /
                 mpzToDD((2-orsa::kronecker(0,m))*(2*l+1)*orsa::factorial(l-m)));
        // stored[l][m] = val;
        return val;
    }
    
public:
    class FiveVars {
    public:
        int tau, l, m, u, nu;
        // T Q;
        T ABQ_R0;
    };
    
    inline static dd_real Q(const FiveVars & fv) {
        return
            mpzToDD(orsa::power_sign(fv.u+fv.nu) *
                    orsa::binomial(fv.m,2*fv.nu+fv.tau) *
                    orsa::binomial(fv.l,fv.u) *
                    orsa::binomial(2*fv.l-2*fv.u,fv.l) *
                    orsa::pochhammer(mpz_class(fv.l-fv.m-2*fv.u+1),fv.m)) 
            / mpzToDD(orsa::int_pow((mpz_class)2,fv.l));
    }

public:
    class Trii {
    public:
        int nx, ny, nz;
        size_t index;
    };
    
protected:
    // check if it is in the SQLite db
    bool inDB(double & val,
              const size_t & index) const {
        char **result;
        int nrows, ncols;
        char * zErr;
        char sql_line[1024];
        sprintf(sql_line,
                "SELECT * FROM simplex WHERE id=%zi",
                index);
        int rc = sqlite3_get_table(db,sql_line,&result,&nrows,&ncols,&zErr);
        //
        if (rc != SQLITE_OK) {
            if (zErr != NULL) {
                fprintf(stderr,"SQL error: %s\n",zErr);
                sqlite3_free(zErr);
            }
        }
        // ORSA_DEBUG("nrows: %i  ncols: %i",nrows, ncols);
        //
        /* for (int i=0; i<nrows; ++i) {
           for (int j=0; j<ncols; ++j) {
           // i=0 is the header
           const int index = (i+1)*ncols+j;
           ORSA_DEBUG("result[%i] = %s",index, result[index]);
           }
           }
        */
        //
        bool haveit=false;
        if (nrows==0) {
            // nothing, but must keep this case!
        } else if (nrows==1) {
            // val[index] = atof(result[9]);
            // needToCompute = false;
            val=atof(result[9]);
            haveit=true;
        } else { // if (nrows>size_H) {
            ORSA_ERROR("database corrupted, only 1 entry per index is admitted");
        }
        //
        sqlite3_free_table(result);
        return haveit;
    }
    
public:
    double getIntegral(const int & nx, const int & ny, const int & nz, const bool verbose=true) const {
        const size_t degree = nx+ny+nz;
        const size_t index = getIndex(nx,ny,nz);
        if (val.size() <= index) {
            val.resize(index+1);
        }
        if (!val[index].isSet()) {

            double DBval;
            const bool needToCompute = !inDB(DBval,index);
            if (!needToCompute) {
                val[index] = DBval;
            } else {
                if (verbose) ORSA_DEBUG("value for [%i][%i][%i] not available, computing it...",nx,ny,nz);
                
                std::vector<FiveVars> fvv;
                
                // size_t num_skipped=0;
                
                const size_t l_max = std::max(norm_A.size()-1,norm_B.size()-1);
                for (size_t tau=0; tau<=1; ++tau) {
                    for (size_t l=0; l<=l_max; ++l) {
                        for (size_t m=0; m<=l; ++m) {
                            for (size_t u=0; u<=(l/2); ++u) {
                                if (m>=tau) {
                                    for (size_t nu=0; nu<=((m-tau)/2); ++nu) {
                                        
                                        FiveVars fv;
                                        fv.tau=tau;
                                        fv.l=l;
                                        fv.m=m;
                                        fv.u=u;
                                        fv.nu=nu;
                                        // fv.Q = Q(fv);
                                        // NOTE the normalization factor included here...
                                        // fv.ABQ_R0 = orsa::int_pow(norm_A[l][m],1-tau)*orsa::int_pow(norm_B[l][m],tau)*Q(fv)/normalization_factor(l,m);
                                        // fv.ABQ_R0 = pow(norm_A[l][m],1-tau)*pow(norm_B[l][m],tau)*Q(fv)/normalization_factor(l,m);
                                        // dividing here by R0 to obtain a global 1/R0^Nr
                                        fv.ABQ_R0 = oneOverR0*pow(norm_A[l][m],1-tau)*pow(norm_B[l][m],tau)*Q(fv)/normalization_factor(l,m);
                                        if (fv.ABQ_R0 != 0.0) {
                                            fvv.push_back(fv);
                                        } else {
                                            // ++num_skipped;
                                            // ORSA_DEBUG("skipped!");
                                        }
                                        
                                        /* ORSA_DEBUG("tau: %i   l: %i   m: %i   u: %i   nu: %i    Q: %g",
                                           tau,l,m,u,nu,
                                           ::to_double(fv.Q));
                                        */
                                    }

                                    
                                }
                            }
                        }
                    }
                }
                // ORSA_DEBUG("fvv.size: %i   skipped: %i",fvv.size(),num_skipped);
                
                const size_t Nr = nx+ny+nz+3;
                
                // include those integrals needed at same Nr
                std::vector<Trii> ii; // includedIntegrals;
                Trii trii;
                trii.nx=nx; trii.ny=ny; trii.nz=nz; trii.index=getIndex(nx,ny,nz);
                ii.push_back(trii);
                for (int gx=0; gx<=degree; ++gx) {
                    for (int gy=0; gy<=degree-gx; ++gy) {
                        // for (int gz=0; gz<=degree-gx-gy; ++gz) {
                        int gz = degree-gx-gy;
                        if (gx==nx && gy==ny && gz==nz) continue; // already included this one...
                        // ORSA_DEBUG("testing: %i %i %i   deg: %i",gx,gy,gz,degree);
                        if ((gx+gy+gz)==degree) {
                            trii.nx=gx;
                            trii.ny=gy;
                            trii.nz=gz;
                            trii.index=getIndex(gx,gy,gz);
                            if (val.size() <= trii.index) {
                                val.resize(trii.index+1);
                            }
                            if ( (!val[trii.index].isSet()) &&
                                 (!inDB(DBval,trii.index)) ) {
                                ii.push_back(trii);
                                if (verbose) ORSA_DEBUG("also including computation for [%i][%i][%i]...",gx,gy,gz);
                            }
                        }
                        // }
                    }
                }
                
                // dd_real big_sum = 0;
                std::vector<dd_real> big_sum;
                big_sum.resize(ii.size());
                for (size_t jj=0; jj<ii.size(); ++jj) {
                    big_sum[jj] = 0.0;
                }
                
                std::vector<size_t> pos;
                pos.resize(Nr);
                
                while (1) {

                    if (0) {
                        // debug only...
                        std::cout << "pos:";
                        size_t p=pos.size();
                        while(p!=0) {
                            --p;
                            std::cout << " " << pos[p];
                        }
                        std::cout << std::endl;
                    }
                    
                    if (0) {
                        // debug only...
                        {
                            size_t p=pos.size();
                            while(p!=0) {
                                --p;
                                const FiveVars & fv = fvv[pos[p]];
                                char line[4096];
                                gmp_sprintf(line," [%i,%i,%i,%i,%i|%g]",fv.tau,fv.l,fv.m,fv.u,fv.nu,::to_double(fv.ABQ_R0));
                                std::cout << line; // << " ";
                                // if (p=!0) cout << " ";
                            }
                        }
                        std::vector<size_t> count;
                        count.resize(fvv.size());
                        for (size_t p=0; p<pos.size(); ++p) {
                            ++count[pos[p]];
                        }
                        mpz_class factor = orsa::factorial(Nr);
                        for (size_t c=0; c<count.size(); ++c) {
                            factor /= orsa::factorial(count[c]);
                        }
                        // total += factor;
                        // cout << "x " << factor.get_mpz_t() << endl;
                        char cn[1024];
                        gmp_sprintf(cn," %Zi",factor.get_mpz_t());
                        std::cout << cn;
                        std::cout << std::endl;
                    }
                    
                    std::vector<int> count;
                    count.resize(fvv.size());
                    for (size_t c=0; c<count.size(); ++c) {
                        count[c] = 0;
                    }
                    for (size_t s=0; s<Nr; ++s) {
                        ++count[pos[s]];
                    }
                    mpz_class binomial_factor = orsa::factorial(Nr);
                    for (size_t c=0; c<count.size(); ++c) {
                        binomial_factor /= orsa::factorial(count[c]);
                    }
                    
                    // ORSA_DEBUG("binomial_factor: %Zi",binomial_factor.get_mpz_t());
                    
                    // dd_real just_norm_factors = 1.0;
                    dd_real coefficients_factor = 1.0;
                    //  dividing here already...
                    // dd_real coefficients_factor = pow(oneOverR0,Nr)/Nr; // 1/(Nr*R0^Nr)
                    for (size_t c=0; c<fvv.size(); ++c) {
                        const FiveVars & fv = fvv[c];
                        /* coefficients_factor *=
                           int_pow(int_pow(norm_A[fv.l][fv.m],1-fv.tau)*
                           int_pow(norm_B[fv.l][fv.m],fv.tau),
                           count[c]) *
                           int_pow(fv.Q,count[c])
                           / int_pow(normalization_factor(fv.l,fv.m),count[c]); // dealing with normalized coefficients...
                        */
                        coefficients_factor *= pow(fv.ABQ_R0,count[c]);
                        // just_norm_factors *= int_pow(normalization_factor(fv.l,fv.m),count[c]);
                        // if (coefficients_factor != 0) ORSA_DEBUG("nf[%i][%i] = %g    cf: %g",fv.l,fv.m,::to_double(normalization_factor(fv.l,fv.m)),::to_double(coefficients_factor));
                    }
                    
                    // ORSA_DEBUG("coefficients_factor: %g",::to_double(coefficients_factor));
                    
                    if (coefficients_factor != 0.0) {
                        
                        // initial values...
                        /* int pow_cos_phi=nx;
                           int pow_sin_phi=ny;
                           int pow_cos_theta=nz;
                           int pow_sin_theta=nx+ny+1;
                        */
                        int pow_cos_phi=0;
                        int pow_sin_phi=0;
                        int pow_cos_theta=0;
                        int pow_sin_theta=0;
                        for (int c=0; c<fvv.size(); ++c) {
                            /* ORSA_DEBUG("%i %i %i %i",
                               pow_cos_phi,
                               pow_sin_phi,
                               pow_cos_theta,
                               pow_sin_theta); */
                            const FiveVars & fv = fvv[c];
                            pow_cos_phi   += count[c]*(fv.m-2*fv.nu-fv.tau);
                            pow_sin_phi   += count[c]*(2*fv.nu+fv.tau);
                            pow_cos_theta += count[c]*(fv.l-fv.m-2*fv.u);
                            pow_sin_theta += count[c]*(fv.m);
                        }

                        for (size_t jj=0; jj<ii.size(); ++jj) {
                        
                            /* const dd_real   factor_phi_integral = integral_csk(pow_cos_phi,pow_sin_phi,2);
                               const dd_real factor_theta_integral = integral_csk(pow_cos_theta,pow_sin_theta,1);
                            */
                            
                            const dd_real   factor_phi_integral = integral_csk(pow_cos_phi+ii[jj].nx,
                                                                               pow_sin_phi+ii[jj].ny,
                                                                               2);
                            
                            const dd_real factor_theta_integral = integral_csk(pow_cos_theta+ii[jj].nz,
                                                                               pow_sin_theta+ii[jj].nx+ii[jj].ny+1,
                                                                               1);
                            
                            // ORSA_DEBUG("csk: %i %i %i  itg: %g",pow_cos_phi,pow_sin_phi,2,::to_double(factor_phi_integral));
                            // ORSA_DEBUG("csk: %i %i %i  itg: %g",pow_cos_theta,pow_sin_theta,1,::to_double(factor_theta_integral));
                            
                            // const T old_big_sum = big_sum;
                            big_sum[jj] +=
                                mpzToDD(binomial_factor) *
                                coefficients_factor *
                                factor_phi_integral *
                                factor_theta_integral;
                            
                            // ORSA_DEBUG("big_sum[%i] = %g",jj,::to_double(big_sum[jj]));
                            
                            /* ORSA_DEBUG("/ff/ %g %g %g %g",
                               binomial_factor.get_d(),
                               ::to_double(coefficients_factor),
                               ::to_double(factor_phi_integral),
                               ::to_double(factor_theta_integral));
                            */
                            
                        }
                        
                        // ORSA_DEBUG("big_sum: %g",::to_double(big_sum));
                        
                        /* 
                           if (0) if (big_sum != old_big_sum) {
                           // debug only...
                           {
                           size_t p=pos.size();
                           while(p!=0) {
                           --p;
                           const FiveVars & fv = fvv[pos[p]];
                           char line[4096];
                           gmp_sprintf(line," [%i,%i,%i,%i,%i|%g]",fv.tau,fv.l,fv.m,fv.u,fv.nu,::to_double(fv.ABQ_R0));
                           std::cout << line; // << " ";
                           // if (p=!0) cout << " ";
                           }
                           }
                           char line[4096];
                           gmp_sprintf(line," s: %g   ds: %g   bin: %g   coeff: %g   phi: %g   theta: %g   jnf: %g",
                           ::to_double(big_sum),
                           ::to_double(big_sum-old_big_sum),                                       
                           binomial_factor.get_d(),
                           ::to_double(coefficients_factor),
                           ::to_double(factor_phi_integral),
                           ::to_double(factor_theta_integral),
                           ::to_double(just_norm_factors));
                           std::cout << line << std::endl;
                           }
                        */
                        
                    }
                    
                    // now, increment while avoiding repetitions
                    bool changed=false;
                    for (size_t p=0; p<pos.size(); ++p) {
                        if (pos[p]<fvv.size()-1) {
                            ++pos[p];
                            changed=true;
                            for (size_t s=0; s<p; ++s) {
                                if (pos[s]>pos[p]) {
                                    pos[s]=pos[p];
                                }
                            }
                            break;
                        }
                    }
                    if (!changed) break; // done
                    
                }
                for (size_t jj=0; jj<ii.size(); ++jj) {
                    big_sum[jj] /= Nr;
                    // big_sum[jj] /= pow(R0,Nr); // dividing earlier, in fv.ABQ_R0
                }
                
                // val[index] = aux_02(retVal,orsa::pochhammer(mpz_class(N+1),degree));
                for (size_t jj=0; jj<ii.size(); ++jj) {
                    // val[index] = ::to_double(big_sum);
                    val[ii[jj].index] = ::to_double(big_sum[jj]);
                }
                
                for (size_t jj=0; jj<ii.size(); ++jj) {
                    
                    /* ORSA_DEBUG("saving on DB element [%i][%i][%i] = %g   index: %i",
                       ii[jj].nx,ii[jj].ny,ii[jj].nz,
                       (*val[ii[jj].index]),
                       ii[jj].index);
                    */
                    
                    char * zErr;
                    char sql_line[1024];
                    sprintf(sql_line,
                            "INSERT INTO simplex VALUES(%zi,%zi,%zi,%zi,%.12e)",
                            ii[jj].index,ii[jj].nx,ii[jj].ny,ii[jj].nz,
                            (*val[ii[jj].index]));
                    int rc;
                    do {
                        rc = sqlite3_exec(db,sql_line,NULL,NULL,&zErr);
                        if (rc==SQLITE_BUSY) {
                            ORSA_DEBUG("database busy, retrying...");
                            usleep(100000);
                        }
                    } while (rc==SQLITE_BUSY);
                    if (rc != SQLITE_OK) {
                        if (zErr != NULL) {
                            fprintf(stderr,"SQL error: %s\n",zErr);
                            sqlite3_free(zErr);
                        }
                    }
                }
                
            }
        }
        return val[index];
    }
public:
    typedef std::vector< std::vector< std::vector<size_t> > > IndexTableType;
    static IndexTableType indexTable;
public:
    static void updateIndexTable(const size_t & requestedDegree) {
        if (indexTable.size() < (requestedDegree+1)) {
            indexTable.resize(requestedDegree+1);
            for (size_t i=0; i<=requestedDegree; ++i) {
                indexTable[i].resize(requestedDegree+1-i);
                for (size_t j=0; j<=requestedDegree-i; ++j) {
                    indexTable[i][j].resize(requestedDegree+1-i-j);
                }
            }
            size_t idx=0;
            size_t degree=0;
            while (degree <= requestedDegree) {
                for (unsigned int i=0; i<=degree; ++i) {
                    for (unsigned int j=0; j<=degree; ++j) {
                        for (unsigned int k=0; k<=degree; ++k) {
                            if (i+j+k==degree) {
                                // ORSA_DEBUG("inserting %i-%i-%i  index: %i",i,j,k,idx);
                                indexTable[i][j][k] = idx++;
                            }
                        }
                    }
                }
                ++degree;
            }
        }
    }
public:
    static size_t getIndex(const size_t & nx, const size_t & ny, const size_t & nz) {
        const size_t requestedDegree=nx+ny+nz;
        updateIndexTable(requestedDegree);
        return indexTable[nx][ny][nz]; 
    }
    
    // same as above, but for 4 indices
public:
    typedef std::vector< std::vector< std::vector< std::vector<size_t> > > > Index4TableType;
    static Index4TableType index4Table;
public:
    static void updateIndex4Table(const size_t & requestedDegree) {
        if (index4Table.size() < (requestedDegree+1)) {
            index4Table.resize(requestedDegree+1);
            for (size_t i=0; i<=requestedDegree; ++i) {
                index4Table[i].resize(requestedDegree+1-i);
                for (size_t j=0; j<=requestedDegree-i; ++j) {
                    index4Table[i][j].resize(requestedDegree+1-i-j);
                    for (size_t k=0; k<=requestedDegree-i-j; ++k) {
                        index4Table[i][j][k].resize(requestedDegree+1-i-j-k);
                    }
                }
            }
            size_t idx=0;
            size_t degree=0;
            while (degree <= requestedDegree) {
                for (unsigned int i=0; i<=degree; ++i) {
                    for (unsigned int j=0; j<=degree; ++j) {
                        for (unsigned int k=0; k<=degree; ++k) {
                            for (unsigned int l=0; l<=degree; ++l) {
                                if (i+j+k+l==degree) {
                                    // ORSA_DEBUG("inserting %i-%i-%i-%i  index4: %i",i,j,k,l,idx);
                                    index4Table[i][j][k][l] = idx++;
                                }
                            }
                        }
                    }
                }
                ++degree;
            }
        }
    }
public:
    static size_t getIndex4(const size_t & i, const size_t & j, const size_t & k, const size_t & l) {
        const size_t requestedDegree=i+j+k+l;
        updateIndex4Table(requestedDegree);
        return index4Table[i][j][k][l]; 
    }
};

#endif // _SH2IJK_H_
