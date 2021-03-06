#ifndef _SH2IJK_H_
#define _SH2IJK_H_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/double.h>
#include <orsa/shape.h>
#include <orsa/unit.h>
#include <orsa/util.h>

// SQLite3
#include "sqlite3.h"

#include <libgen.h>

// #include <gsl/gsl_sf_gamma.h>

#include "simplex.h"
#include "mpreal.h"

inline std::string getSqliteDBFileName_SH(const std::string & inputFile,
                                          const double & R0) {
    char line[4096];
    char arg_of_dirname[4096];
    char arg_of_basename[4096];
    strcpy(arg_of_dirname,inputFile.c_str());
    strcpy(arg_of_basename,inputFile.c_str());
    sprintf(line,"%s/.%s_simplex_%gkm.sqlite",
            dirname(arg_of_dirname),
            basename(arg_of_basename),
            orsa::FromUnits(R0,orsa::Unit::KM,-1));
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

typedef std::vector< std::vector<double> > SHcoeff;

inline void writeSH(const SHcoeff & norm_A,
                    const SHcoeff & norm_B,
                    const std::string & outputFile,
                    const orsa::Unit::LengthUnit & outputFileLengthScale = orsa::Unit::METER) {
    
    FILE * fp = fopen(outputFile.c_str(),"w");
    
    if (!fp) {
        ORSA_DEBUG("cannot open file [%s]",outputFile.c_str());
        return;
    }
    
    for (size_t l=0; l<norm_A.size(); ++l) {
        for (size_t m=0; m<=l; ++m) {
            if (m!=0) {
                gmp_fprintf(fp,"%3i %3i %+16.6f %+16.6f\n",
                            l,
                            m,
                            orsa::FromUnits(norm_A[l][m],outputFileLengthScale,-1),
                            orsa::FromUnits(norm_B[l][m],outputFileLengthScale,-1));
            } else {
                gmp_fprintf(fp,"%3i %3i %+16.6f\n",
                            l,
                            m,
                            orsa::FromUnits(norm_A[l][m],outputFileLengthScale,-1));
            }
        }
    }
    
    fclose(fp);
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
    char s_l[4096], s_m[4096], s_A[4096], s_B[4096];
    while (fgets(line,4096,fp)) {
        int num_good = sscanf(line,
                              "%s %s %s %s",
                              s_l, s_m, s_A, s_B);
        // ORSA_DEBUG("[%s] [%s] [%s] [%s] (num_good=%i)",s_l, s_m, s_A, s_B, num_good);
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
            
            norm_A[l][m] = orsa::FromUnits(atof(s_A),inputFileLengthScale);
            if ( (num_good==4) && (m>0) ) {
                norm_B[l][m] = orsa::FromUnits(atof(s_B),inputFileLengthScale);
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

template <typename F> class SHIntegration : public osg::Referenced {
public:
    SHIntegration(const SHcoeff & nA,
                  const SHcoeff & nB,
                  const F & R0_,
                  // const F & epsabs_,
                  // const F & epsrel_,
                  const std::string & SQLiteDBFileName) :
        osg::Referenced(),
        norm_A(nA),
        norm_B(nB),
        oneOverR0(1.0/R0_) {
        // epsabs(epsabs_),
        // epsrel(epsrel_) {
        
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
    /*
    void reserve(const size_t maxDegree) {
        const size_t maxIndex = std::max(getIndex(maxDegree,0,0),
                                         std::max(getIndex(0,maxDegree,0),
                                                  getIndex(0,0,maxDegree)));
        val.reserve(maxIndex+1);
    }
    */
protected:
    const SHcoeff & norm_A;
    const SHcoeff & norm_B;
    const F oneOverR0;
    // const F epsabs;
    // const F epsrel;
    mutable std::vector< orsa::Cache<double> > val; // integral value
protected:
    // needed to work with SQLite database
    sqlite3 * db;
protected:
    // n!!
    static mpfr::mpreal bifactorial(int n) {
        if (n<2) {
            return 1;
        } else {
            mpfr::mpreal retVal=n;
            while (n>3) {
                n -= 2;
                retVal *= n;
            }
            return retVal;
        }
    }
    // demiGamma(n/2)
    static mpfr::mpreal demiGamma(const int n) {
        static const mpfr::mpreal sqrt_pi = sqrt(mpfr::const_pi());
        static const mpfr::mpreal one_over_sqrt_2 = mpfr::mpreal(1)/sqrt(mpfr::mpreal(2));
        return sqrt_pi*bifactorial(n-2)*pow(one_over_sqrt_2,n-1);
    }
    // fast!
    /*
    static mpfr::mpreal triple_gamma(const int c, const int s) {
        return (demiGamma(c+1)*demiGamma(s+1)/demiGamma(c+s+2));
    }
    */
    // fast + lookup table
    static mpfr::mpreal triple_gamma(const unsigned int c, const unsigned int s) {
        static std::vector<std::vector<orsa::Cache<mpfr::mpreal> > > table;
        if (table.size()>c) {
            if (table[c].size()>s) {
                if (table[c][s].isSet()) {
                    return table[c][s];
                } else {
                    table[c][s] = (demiGamma(c+1)*demiGamma(s+1)/demiGamma(c+s+2));
                    return table[c][s];
                }
            } else {
                table[c].resize(s+1);
                return triple_gamma(c,s);
            }
        } else {
            table.resize(c+1);
            return triple_gamma(c,s);
        }
    }
    /*
    // slow!
    mpfr::mpreal triple_gamma(const mpfr::mpreal c, const mpfr::mpreal s) const {
        return (mpfr::gamma((c+1)/2)*mpfr::gamma((s+1)/2)/mpfr::gamma((c+s+2)/2));
    }
    */
protected:
    // 5-10% faster than version without lookup table...
    mpfr::mpreal integral_csk(const int & c, const int & s, const int & k) const {
        static std::vector<std::vector<std::vector<orsa::Cache<mpfr::mpreal> > > > table;
        if ((c<0) || (s<0) || (k<0)) return 0; // tested only for natural numbers values
        if (table.size()>c) {
            if (table[c].size()>s) {
                if (table[c][s].size()>k) {
                    if (table[c][s][k].isSet()) {
                        return table[c][s][k];
                    } else {
                        // calculating it here...
                        const bool c_even = (c%2==0);
                        if (c_even) {
                            const bool s_even = (s%2==0);
                            if (s_even) {
                                table[c][s][k] = k*triple_gamma(c,s);
                                return table[c][s][k];
                            } else {
                                const bool k_even = (k%2==0);
                                if (k_even) {
                                    table[c][s][k] = 0;
                                    return table[c][s][k];
                                } else {
                                    table[c][s][k] = triple_gamma(c,s);
                                    return table[c][s][k];
                                }
                            }
                        } else {
                            table[c][s][k] = 0;
                            return table[c][s][k];
                        }
                    }
                } else {
                    table[c][s].resize(k+1);
                    return integral_csk(c,s,k);
                }
            } else {
                table[c].resize(s+1);
                return integral_csk(c,s,k);
            }
        } else {
            table.resize(c+1);
            return integral_csk(c,s,k);
        }
    }
    /*
    mpfr::mpreal integral_csk(const int & c, const int & s, const int & k) const {
        if ((c<0) || (s<0) || (k<0)) return 0; // tested only for natural numbers values
        const bool c_even = (c%2==0);
        if (c_even) {
            const bool s_even = (s%2==0);
            if (s_even) {
                return k*triple_gamma(c,s);
            } else {
                const bool k_even = (k%2==0);
                if (k_even) {
                    return 0;
                } else {
                    return triple_gamma(c,s);
                }
            }
        } else {
            return 0;
        }
    }
    */
protected:
    mpfr::mpreal normalization_factor(const unsigned long int l, const unsigned long int m) const {
        return sqrt(mpfr::fac_ui(l+m)/((2-kronecker(0,m))*(2*l+1)*mpfr::fac_ui(l-m)));
    }
    
public:
    class FiveVars {
    public:
        int tau, l, m, u, nu;
        // T Q;
        F ABQ_R0;
    public:
        // note: sorting by absolute value!
        static bool sort_by_absolute_larger_to_smaller(const FiveVars & x, const FiveVars & y) {
            return (fabs(x.ABQ_R0) > fabs(y.ABQ_R0));
        }
    };
    inline static F Q(const FiveVars & fv) {                    
        return (orsa::power_sign(fv.u+fv.nu)*binomial(fv.m,2*fv.nu+fv.tau)*binomial(fv.l,fv.u)*binomial(2*fv.l-2*fv.u,fv.l)*pochhammer(fv.l-fv.m-2*fv.u+1,fv.m))/(mpfr::pow(mpfr::mpreal(2),fv.l));
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
                
                // ORSA_DEBUG("NOTE: using accuracy coefficients epsabs = %g and epsrel = %g",::to_double(epsabs),::to_double(epsrel));
                
                const size_t l_max = std::max(norm_A.size()-1,norm_B.size()-1);
                
                std::vector<FiveVars> fvv;
                // fvv.reserve(l_max*l_max*l_max);
                
                for (size_t tau=0; tau<=1; ++tau) {
                    for (size_t l=0; l<=l_max; ++l) {
                        for (size_t m=0; m<=l; ++m) {
                            const F nf = normalization_factor(l,m);
                            for (size_t u=0; u<=(l/2); ++u) {
                                if (m>=tau) {
                                    for (size_t nu=0; nu<=((m-tau)/2); ++nu) {
                                        
                                        FiveVars fv;
                                        fv.tau=tau;
                                        fv.l=l;
                                        fv.m=m;
                                        fv.u=u;
                                        fv.nu=nu;
                                        // NOTE the normalization factor included here...
                                        fv.ABQ_R0 = oneOverR0*pow(norm_A[l][m],1.0-tau)*pow(norm_B[l][m],tau)*Q(fv)/nf;
                                        // ORSA_DEBUG("%i %i %i %i %i Q: %g   nf: %g",tau,l,m,u,nu,Q(fv).toDouble(),normalization_factor(l,m).toDouble());
                                        if (fv.ABQ_R0 != 0) {
                                            fvv.push_back(fv);
                                            // if (verbose) ORSA_DEBUG("[%i,%i,%i,%i,%i] = %+12.9f",fv.tau,fv.l,fv.m,fv.u,fv.nu,::to_double(fv.ABQ_R0));
                                        } 
                                    }
                                }
                            }
                        }
                    }
                }
                
                std::sort(fvv.begin(),fvv.end(),FiveVars::sort_by_absolute_larger_to_smaller);
                if (verbose) for (size_t c=0; c<fvv.size(); ++c) {
                    const FiveVars & fv = fvv[c];
                    ORSA_DEBUG("[%i,%i,%i,%i,%i] = %+12.9f",fv.tau,fv.l,fv.m,fv.u,fv.nu,fv.ABQ_R0.toDouble());
                }
                
                ORSA_DEBUG("fvv.size(): %u",fvv.size());
                
                const size_t Nr = nx+ny+nz+3;
                
                // include those integrals needed at same Nr
                std::vector<Trii> ii; // includedIntegrals;
                Trii trii;
                trii.nx=nx; trii.ny=ny; trii.nz=nz; trii.index=getIndex(nx,ny,nz);
                ii.push_back(trii);
                for (int gx=0; gx<=degree; ++gx) {
                    for (int gy=0; gy<=degree-gx; ++gy) {
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
                    }
                }
                
                std::vector<F> big_sum;
                big_sum.resize(ii.size());
                for (size_t jj=0; jj<ii.size(); ++jj) {
                    big_sum[jj] = 0;
                }
                
                std::vector<size_t> pos;
                pos.resize(Nr);
                for (size_t p=0; p<Nr; ++p) {
                    pos[p] = 0;
                }
                
                // orsa::Cache<F> min_abs_big_sum;
                // orsa::Cache<F> coefficients_factor_threshold;
                
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
                                gmp_sprintf(line," [%i,%i,%i,%i,%i|%g]",fv.tau,fv.l,fv.m,fv.u,fv.nu,fv.ABQ_R0.toDouble());
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
                    
                    // bool skip_term=false;
                    
                    F binomial_factor = mpfr::fac_ui(Nr);
                    for (size_t c=0; c<count.size(); ++c) {
                        binomial_factor /= mpfr::fac_ui(count[c]);
                    }
                    
                    F coefficients_factor = binomial_factor;
                    for (size_t c=0; c<fvv.size(); ++c) {
                        coefficients_factor *= pow(fvv[c].ABQ_R0,count[c]);
                        
                        /*
                        // can skip this term?
                        if ((nx==0) && (ny==0) && (nz==0)) {
                            // cannot skip for 0,0,0 term, it would introduce a larger error
                        } else if (coefficients_factor_threshold.isSet()) {
                            if ( (fabs(coefficients_factor) < coefficients_factor_threshold) && (fabs(fvv[c].ABQ_R0) < 1.0) ) {
                                
                                if (verbose) {
                                    if (1) {
                                        std::cout << "pos:";
                                        size_t p=pos.size();
                                        while(p!=0) {
                                            --p;
                                            std::cout << " " << pos[p];
                                        }
                                        std::cout << std::endl;
                                    }
                                    ORSA_DEBUG("skipping term, min_abs_big_sum = %g   coefficients_factor: %g",(*min_abs_big_sum).toDouble(),coefficients_factor.toDouble());
                                }
                                
                                skip_term=true;
                                break;
                            }
                        }
                        */
                    }
                    
                    // if (verbose) ORSA_DEBUG("coefficients_factor: %g",coefficients_factor.toDouble());
                    
                    /* 
                    if ( (!skip_term) &&
                         (coefficients_factor != 0) ) {
                    */
                    if (coefficients_factor != 0) {
                            
                        // initial values...
                        int pow_cos_phi=0;
                        int pow_sin_phi=0;
                        int pow_cos_theta=0;
                        int pow_sin_theta=0;
                        for (int c=0; c<fvv.size(); ++c) {
                            const FiveVars & fv = fvv[c];
                            pow_cos_phi   += count[c]*(fv.m-2*fv.nu-fv.tau);
                            pow_sin_phi   += count[c]*(2*fv.nu+fv.tau);
                            pow_cos_theta += count[c]*(fv.l-fv.m-2*fv.u);
                            pow_sin_theta += count[c]*(fv.m);
                            
                        /*  ORSA_DEBUG("%i %i %i %i",
                            pow_cos_phi,
                            pow_sin_phi,
                            pow_cos_theta,
                            pow_sin_theta); */
                            
                        }

                        for (size_t jj=0; jj<ii.size(); ++jj) {
                            const F factor_phi_integral   = integral_csk(pow_cos_phi+ii[jj].nx,pow_sin_phi+ii[jj].ny,2);
                            const F factor_theta_integral = integral_csk(pow_cos_theta+ii[jj].nz,pow_sin_theta+ii[jj].nx+ii[jj].ny+1,1);
                            
                            /* ORSA_DEBUG("jj: %02i   csk: %i %i %i   itg: %g",
                               jj,
                               pow_cos_phi+ii[jj].nx,
                               pow_sin_phi+ii[jj].ny,
                               2,
                               ::to_double(factor_phi_integral));
                               ORSA_DEBUG("jj: %02i   csk: %i %i %i   itg: %g",
                               jj,
                               pow_cos_theta+ii[jj].nz,
                               pow_sin_theta+ii[jj].nx+ii[jj].ny+1,
                               1,
                               ::to_double(factor_theta_integral));
                            */
                            
                            big_sum[jj] +=
                                coefficients_factor *
                                factor_phi_integral *
                                factor_theta_integral;
                            
                            // ORSA_DEBUG("big_sum[%i] = %40.20f",jj,big_sum[jj].toDouble());

                            /* 
                               if (verbose) ORSA_DEBUG("/ff/ %6g %+16.6f %6.3f %6.3f    big_sum[%i] = %+12.9g",
                               binomial_factor.get_d(),
                               ::to_double(coefficients_factor),
                               ::to_double(factor_phi_integral),
                               ::to_double(factor_theta_integral),
                               jj,::to_double(big_sum[jj]));
                            */
                            
                        }
                        
                        /*
                        // important update
                        min_abs_big_sum.reset();
                        for (size_t jj=0; jj<ii.size(); ++jj) {
                            if (big_sum[jj] != 0.0) {
                                min_abs_big_sum.setIfSmaller(fabs(big_sum[jj]));
                            }
                        }
                        // the factor (4*pi) comes from the two integrals which are up to 2*pi and 2, respectively
#warning check again the maximum value of each integral!
                        if (min_abs_big_sum.isSet()) {
                            // coefficients_factor_threshold = epsrel*min_abs_big_sum/(4*orsa::pi()*Nr_factorial_dd);
                            coefficients_factor_threshold = (epsabs + epsrel*min_abs_big_sum) / (4*orsa::pi());
                        } else {
                            coefficients_factor_threshold.reset();
                        }
                        */
                    }
                    
                    // now, increment while avoiding repetitions
                    /*
                    if (skip_term) {
                        bool changed=false;
                        for (size_t p=0; p<pos.size(); ++p) {
                            if (pos[p]<fvv.size()-1) {
                                if (p>0) {
                                    if (pos[p]!=pos[p-1]) {
                                        ++pos[p];
                                        changed=true;
                                        for (size_t s=0; s<p; ++s) {
                                            pos[s]=pos[p];
                                        }
                                        break;
                                    }
                                }
                                if (p<pos.size()-1) {
                                    if (pos[p]!=pos[p+1]) {
                                        ++pos[p+1];
                                        changed=true;
                                        for (size_t s=0; s<=p; ++s) {
                                            pos[s]=pos[p+1];
                                        }
                                        break;
                                    }
                                }
                            }
                        }
                        if (!changed) break; // done
                    } else {
                    */
                    
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
                    
                    // }
                    
                }
                for (size_t jj=0; jj<ii.size(); ++jj) {
                    big_sum[jj] /= Nr;
                }
                
                for (size_t jj=0; jj<ii.size(); ++jj) {
                    val[ii[jj].index] = big_sum[jj].toDouble();
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
