#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/double.h>
#include <orsa/shape.h>
#include <orsa/unit.h>

#include <qd/dd_real.h>
#include <qd/qd_real.h>

#warning if shape is strongly concave and a simplex covers volume outside the body shape, then the results are incorrect (including volume computations...)

#warning default origin for 4th simplex vertex, should be a parameter of the class??                

// SQLite3
#include "sqlite3.h"

#include <gsl/gsl_sf_gamma.h>

inline double integral_csk_util(const size_t & c, const size_t & s) {
    return gsl_sf_gamma(0.5*(c+1))*gsl_sf_gamma(0.5*(s+1))/gsl_sf_gamma(0.5*(c+s+2));
}
// integral between 0 and k*pi of cos^c(x) sin^s(x) dx
inline double integral_csk(const size_t & c, const size_t & s, const size_t & k) {
    const bool c_even = (c%2==0);
    if (c_even) {
        const bool s_even = (s%2==0);
        if (s_even) {
            return k*integral_csk_util(c,s);
        } else {
            const bool k_even = (k%2==0);
            if (k_even) {
                return 0.0;
            } else {
                return integral_csk_util(c,s);
            }
        }
    } else {
        return 0.0;
    }
}

#warning should change the _simplex_ part of the name to _ijk_ or something like that
inline std::string getSqliteDBFileName(const std::string & inputFile,
                                       const double & R0) {
    char line[1024];
    sprintf(line,"%s_simplex_%gkm.sqlite",inputFile.c_str(),orsa::FromUnits(R0,orsa::Unit::KM,-1));
    return line;
}

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
        ORSA_DEBUG("[%s] [%s] [%s] [%s] (num_good=%i)",s_l, s_m, s_C, s_S, num_good);
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
    
    fclose(fp);
}

template <typename T> class SHIntegration : public osg::Referenced {
public:
    SHIntegration(const SHcoeff & nA,
                  const SHcoeff & nB,
                  const double & R0_,
                  const std::string & SQLiteDBFileName) :
        osg::Referenced(),
        // N(3),
        norm_A(nA),
        norm_B(nB),
        // triShape(s),
        R0(R0_),
        oneOverR0(1.0/R0) {
        
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
    const double R0;
    const double oneOverR0;
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
    /* template <class QD> QD mpzToQD(const mpz_class & z) const {
       char * str = (char *)malloc(mpz_sizeinbase(z.get_mpz_t(),10)+2);
       mpz_get_str(str,10,z.get_mpz_t());
       // ORSA_DEBUG("z: [%Zi]  STR: [%s]",z.get_mpz_t(),str);
       QD x(str);
       free(str);
       return x;
       }
       dd_real mpzToQD(const mpz_class & z) const;
       qd_real mpzToQD(const mpz_class & z) const;
    */
    dd_real mpzToDD(const mpz_class & z) const {
        char * str = (char *)malloc(mpz_sizeinbase(z.get_mpz_t(),10)+2);
        mpz_get_str(str,10,z.get_mpz_t());
        // ORSA_DEBUG("z: [%Zi]  STR: [%s]",z.get_mpz_t(),str);
        dd_real x(str);
        free(str);
        return x;
    }
    qd_real mpzToQD(const mpz_class & z) const {
        char * str = (char *)malloc(mpz_sizeinbase(z.get_mpz_t(),10)+2);
        mpz_get_str(str,10,z.get_mpz_t());
        // ORSA_DEBUG("z: [%Zi]  STR: [%s]",z.get_mpz_t(),str);
        qd_real x(str);
        free(str);
        return x;
    }
    //
    /* template <class U> T to_T(const U & x) const;
       double to_T(const double & x) const { return x; }
       double to_T(const mpf_class & x) const { return x.get_d(); }
       double to_T(const mpz_class & n) const { return n.get_d(); }
    */
    //
    /* double to_double(const double & x) const { return x; }
       double to_double(const mpf_class & x) const { return x.get_d(); }
    */
    //
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
    
    class FiveVars {
    public:
        int tau, l, m, u, nu;
        mpz_class Q;
    };
    
#warning NEED TO TEST THE Q value...
    
    inline static mpz_class Q(const FiveVars & fv) {
        return
            orsa::power_sign(fv.u+fv.nu) *
            orsa::binomial(fv.m,2*fv.nu+fv.tau) *
            orsa::binomial(fv.l,fv.u) *
            orsa::binomial(2*fv.l-2*fv.u,fv.l) *
            orsa::pochhammer(mpz_class(fv.l-fv.m-2*fv.u+1),fv.m) 
            / orsa::int_pow(2,fv.l);
    }
    
public:
    double getIntegral(const size_t & nx, const size_t & ny, const size_t & nz) const {
        const size_t degree = nx+ny+nz;
        const size_t index = getIndex(nx,ny,nz);
        if (val.size() <= index) {
            val.resize(index+1);
        }
        if (!val[index].isSet()) {
            bool needToCompute=true;
            {
                // first check if it is in the SQLite db
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
                if (nrows==0) {
                    // nothing, but must keep this case!
                } else if (nrows==1) {
                    val[index] = atof(result[9]);
                    needToCompute = false;
                } else { // if (nrows>size_H) {
                    ORSA_ERROR("database corrupted, only 1 entry per index is admitted");
                }
                //
                sqlite3_free_table(result);
            }
            
            if (needToCompute) {
                ORSA_DEBUG("value for [%i][%i][%i] not available, computing it...",nx,ny,nz);
                
                std::vector<FiveVars> fvv;
                
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
                                        fv.Q = Q(fv);
                                        fvv.push_back(fv);
                                        
                                        ORSA_DEBUG("tau: %i   l: %i   m: %i   u: %i   nu: %i    Q: %Zi",
                                                   tau,l,m,u,nu,
                                                   fv.Q.get_mpz_t());
                                        
                                    }
                                }
                            }
                        }
                    }
                }
                ORSA_DEBUG("fvv.size: %i",fvv.size());
                
                const size_t Nr = nx+ny+nz+3;
                
                double big_sum = 0.0;
                
                std::vector<size_t> pos;
                pos.resize(Nr);
                
                // double factor = 1.0;
                
                while (1) {
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    {
                        {
                            size_t p=pos.size();
                            while(p!=0) {
                                --p;
                                const FiveVars & fv = fvv[pos[p]];
                                char line[4096];
                                gmp_sprintf(line," [%i,%i,%i,%i,%i|%Zi]",fv.tau,fv.l,fv.m,fv.u,fv.nu,fv.Q.get_mpz_t());
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
                    















                    
                    std::vector<size_t> count;
                    count.resize(fvv.size());
                    for (size_t s=0; s<Nr; ++s) {
                        ++count[pos[s]];
                    }
                    mpz_class binomial_factor = orsa::factorial(Nr);
                    for (size_t c=0; c<count.size(); ++c) {
                        binomial_factor /= orsa::factorial(count[c]);
                    }
                    
                    // factor = binomial_factor.get_d();

                    double coefficients_factor = 1.0;
                    for (size_t c=0; c<fvv.size(); ++c) {
                        const FiveVars & fv = fvv[c];
                        coefficients_factor *=
                            int_pow(int_pow(norm_A[fv.l][fv.m],1-fv.tau)*
                                    int_pow(norm_B[fv.l][fv.m],fv.tau)*
                                    fv.Q.get_d(),count[c]);
                    }
                    
                    // initial values...
                    int pow_cos_phi=nx;
                    int pow_sin_phi=ny;
                    int pow_cos_theta=nz;
                    int pow_sin_theta=nx+ny+1;
                    for (size_t c=0; c<fvv.size(); ++c) {
                        const FiveVars & fv = fvv[c];
                        pow_cos_phi   += count[c]*(fv.m-2*fv.nu-fv.tau);
                        pow_sin_phi   += count[c]*(2*fv.nu-fv.tau);
                        pow_cos_theta += count[c]*(fv.l-fv.m-2*fv.u);
                        pow_sin_theta += count[c]*(fv.m);
                    }
                    const double   factor_phi_integral = integral_csk(pow_cos_phi,pow_sin_phi,2);
                    const double factor_theta_integral = integral_csk(pow_cos_theta,pow_sin_theta,1);
                    
                    big_sum +=
                        binomial_factor.get_d() *
                        coefficients_factor *
                        factor_phi_integral *
                        factor_theta_integral;
                    
                    ORSA_DEBUG("%g %g %g %g",
                               binomial_factor.get_d(),
                               coefficients_factor,
                               factor_phi_integral,
                               factor_theta_integral);
                    
                    ORSA_DEBUG("big_sum: %g",big_sum);
                    
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
                big_sum /= Nr;
                
                // val[index] = aux_02(retVal,orsa::pochhammer(mpz_class(N+1),degree));
                val[index] = big_sum;
                
                char * zErr;
                char sql_line[1024];
                sprintf(sql_line,
                        "INSERT INTO simplex VALUES(%zi,%zi,%zi,%zi,%.12e)",
                        index,nx,ny,nz,
                        (*val[index]));
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

#endif // _SIMPLEX_H_
