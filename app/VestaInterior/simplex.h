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

#include <libgen.h>

inline std::string getSqliteDBFileName_simplex(const std::string & inputFile,
                                               const double & R0) {
    ORSA_DEBUG("test: [%s]",inputFile.c_str());
    char line[4096];
    char arg_of_dirname[4096];
    char arg_of_basename[4096];
    strcpy(arg_of_dirname,inputFile.c_str());
    strcpy(arg_of_basename,inputFile.c_str());
    sprintf(line,"%s/.%s_simplex_%gkm.sqlite",
            dirname(arg_of_dirname),
            basename(arg_of_basename),
            orsa::FromUnits(R0,orsa::Unit::KM,-1));
    ORSA_DEBUG("test: [%s]",inputFile.c_str());
    return line;
}

template <typename T> class SimplexIntegration : public osg::Referenced {
public:
    SimplexIntegration(const orsa::TriShape * s,
                       const double & R0_,
                       const std::string & SQLiteDBFileName) :
        osg::Referenced(),
        N(3),
        triShape(s),
        R0(R0_),
        oneOverR0(1.0/R0) {
        const orsa::TriShape::VertexVector & vv = triShape->getVertexVector();
        const orsa::TriShape::FaceVector   & fv = triShape->getFaceVector();
        aux.resize(fv.size());
        for (size_t fi=0; fi<fv.size(); ++fi) {
            aux[fi].simplexVertexVector.resize(N+1);
            aux[fi].simplexVertexVector[0] = oneOverR0*orsa::Vector(0,0,0);
            aux[fi].simplexVertexVector[1] = oneOverR0*vv[fv[fi].i()];
            aux[fi].simplexVertexVector[2] = oneOverR0*vv[fv[fi].j()];
            aux[fi].simplexVertexVector[3] = oneOverR0*vv[fv[fi].k()];
            //
            // test: is R0 large enough? later we need vars to be between -1 and 1, because they are arguments
            //       of ChebyshevT polynomials, which are orthogonal only in this interval
            if ( (aux[fi].simplexVertexVector[0].length() > 1.0) ||
                 (aux[fi].simplexVertexVector[1].length() > 1.0) ||
                 (aux[fi].simplexVertexVector[2].length() > 1.0) ||
                 (aux[fi].simplexVertexVector[3].length() > 1.0) ) {
                ORSA_DEBUG("problem: R0 value too small");
                exit(0);
            }
            //
            // volume of simplex with the face as base and the origin as 4th vertex
            // if moving 4th point away from origin, more terms must be added
            // aux[fi].volume = fabs((vv[fv[fi].i()]*orsa::externalProduct(vv[fv[fi].j()],vv[fv[fi].k()])) / 6);
            //
            // generic, just in case
            aux[fi].volume =
                fabs((aux[fi].simplexVertexVector[1]-aux[fi].simplexVertexVector[0]) *
                     orsa::externalProduct(aux[fi].simplexVertexVector[2]-aux[fi].simplexVertexVector[0],
                                           aux[fi].simplexVertexVector[3]-aux[fi].simplexVertexVector[0]) / 6.0);
            // ORSA_DEBUG("fi: %02i  volume: %g",fi,(*aux[fi].volume));
        }
        
        int rc = sqlite3_open(SQLiteDBFileName.c_str(),&db);
        /* int rc = sqlite3_open_v2(SQLiteDBFileName.c_str(),
           &db,
           SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,
           "unix-dotfile");
        */
        //
        if (rc) {
            fprintf(stderr,"Can't open db: %s\n",sqlite3_errmsg(db));
            sqlite3_close(db);
            // call exit?
        }
        
        {
            char * zErr;
            // create results table
            std::string sql = "CREATE TABLE if not exists simplex(id INTEGER PRIMARY KEY UNIQUE, nx INTEGER, ny INTEGER, nz INTEGER, integral REAL)";
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
protected:
    virtual ~SimplexIntegration() {
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
    class SimplexInternals {
    public:
        std::vector<orsa::Vector> simplexVertexVector;
        orsa::Cache<T> volume;
    };
protected:
    const size_t N; // N = 3 = dimension of space = number of vertexes in n-dim simplex + 1
    osg::ref_ptr<const orsa::TriShape> triShape;
    const double R0;
    const double oneOverR0;
    mutable std::vector< orsa::Cache<double> > val; // integral value
    mutable std::vector< SimplexInternals > aux;
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
                std::vector<T> val_vol_sum_fun;
                val_vol_sum_fun.resize(degree+1);
                const size_t q_min = (degree==0) ? 0 : 1;
                const orsa::TriShape::FaceVector & fv = triShape->getFaceVector();
                std::vector<size_t> indexVector;
                for (size_t q=q_min; q<=degree; ++q) {
                    T sum_vol_fun = 0;
                    indexVector.resize(q);
                    for (size_t fi=0; fi<fv.size(); ++fi) {
                        for (size_t i=0; i<q; ++i) {
                            indexVector[i] = 0;
                        }
                        T sum_vol_fi = 0;
                        while (1) {
                            // orsa::Vector vI(0,0,0);
                            T vIx = 0;
                            T vIy = 0;
                            T vIz = 0;
                            for (size_t i=0; i<q; ++i) {
                                // vI += aux[fi].simplexVertexVector[indexVector[i]];
                                // ORSA_DEBUG("iV[%02i] = %i",i,indexVector[i]);
                                vIx += aux[fi].simplexVertexVector[indexVector[i]].getX();
                                vIy += aux[fi].simplexVertexVector[indexVector[i]].getY();
                                vIz += aux[fi].simplexVertexVector[indexVector[i]].getZ();
                            }
                            /* sum_vol_fi +=
                               orsa::int_pow(vI.getX(),nx)*
                               orsa::int_pow(vI.getY(),ny)*
                               orsa::int_pow(vI.getZ(),nz);
                            */
                            T vI_pow = 1;
                            for (size_t p=0; p<nx; ++p) { vI_pow *= vIx; }
                            for (size_t p=0; p<ny; ++p) { vI_pow *= vIy; }
                            for (size_t p=0; p<nz; ++p) { vI_pow *= vIz; }
                            sum_vol_fi += vI_pow;
                            /* ORSA_DEBUG("vI = %+g %+g %+g sum term: %+16.6e",vI.getX(),vI.getY(),vI.getZ(),
                               orsa::int_pow(vI.getX(),nx)*
                               orsa::int_pow(vI.getY(),ny)*
                               orsa::int_pow(vI.getZ(),nz));
                            */
                            bool increased = false;
                            for (size_t i=0; i<q; ++i) {
                                if (indexVector[i]<N) {
                                    ++indexVector[i];
                                    increased = true;
                                    for (size_t s=0; s<i; ++s) {
                                        // #warning IMPORTANT: which rule is correct?
                                        indexVector[s] = indexVector[i]; // this one avoids repetitions, i.e. {1,0}, {0,1}
                                        // indexVector[s] = 0; // this one includes repetitions
                                    }
                                    break;
                                }
                            }
                            // ORSA_DEBUG("increased: %i",increased);
                            if (!increased) {
                                break;
                            }                                
                        }
                        sum_vol_fun += sum_vol_fi*(*aux[fi].volume);
                        /* ORSA_DEBUG("degree: %02i  q: %02i sum_vol_fi: %+16.6e  partial sum_vol_fun: %+16.6e",
                           degree,q,sum_vol_fi,sum_vol_fun);
                        */
                    }
                    val_vol_sum_fun[q] = sum_vol_fun;
                }
                T retVal = 0;
                for (size_t q=q_min; q<=degree; ++q) {
                    /* retVal += to_T(mpf_class(orsa::power_sign(degree-q) *
                       orsa::binomial(N+degree,N+q) *
                       val_vol_sum_fun[q]));
                    */
                    /* retVal += to_T(orsa::power_sign(degree-q) *
                       orsa::binomial(N+degree,N+q) *
                       val_vol_sum_fun[q]);
                    */
                    retVal += aux_01(orsa::power_sign(degree-q),
                                     orsa::binomial(N+degree,N+q),
                                     val_vol_sum_fun[q]);
                    /* ORSA_DEBUG("degree: %i  q: %i  factor = binomial(%i,%i): : %g  sign: %i  term: %g",
                       degree,
                       q,
                       N+degree,
                       N+q,
                       sum_vol_fun_factor,sign,
                       val_vol_sum_fun[q]);
                    */
                }
                // val[index] = T(retVal / (orsa::binomial(N+degree,degree)*orsa::factorial(degree))).get_d();
                //
                // the above binomial x factorial can be expressed more compactly with pochhammer(N+1,degree) = pochhammer_Np1[degree]
                /* if (pochhammer_Np1.size() < (degree+1)) {
                   pochhammer_Np1.resize(degree+1);
                   }
                   if (!pochhammer_Np1[degree].isSet()) {
                   pochhammer_Np1[degree] = orsa::pochhammer(mpz_class(N+1),degree);
                   }
                   val[index] = T(retVal / (*pochhammer_Np1[degree])).get_d();
                */
                //
                // or even more simply (slightly slower...)
                // val[index] = to_double((T)(retVal) / (T)(mpf_class(orsa::pochhammer(mpz_class(N+1),degree)).get_d()));
                // val[index] = to_double(retVal / mpf_class(orsa::pochhammer(mpz_class(N+1),degree)));
                val[index] = aux_02(retVal,orsa::pochhammer(mpz_class(N+1),degree));
                
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
