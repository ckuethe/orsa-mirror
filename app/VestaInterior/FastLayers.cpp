#include "FastLayers.h"

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
// typedef dd_real simplex_T;
// typedef qd_real simplex_T;

#warning how to write this using the typedef inside the class?
// template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
// template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;

/*******/

// when using this index, IzzMR2 is used instead of one of the natural Clm or Slm in the gravityData file;
orsa::Cache<size_t> IzzMR2_index;
size_t get_IzzMR2_index() {
    if (IzzMR2_index.isSet()) {
        return IzzMR2_index;
    } else {
        ORSA_DEBUG("error: requesting unset value of IzzMR2_index");
        return -1;
    }
}
orsa::Cache<double> IzzMR2;
double get_IzzMR2() {
    if (IzzMR2.isSet()) {
        return IzzMR2;
    } else {
        ORSA_DEBUG("error: requesting unset value of IzzMR2");
        return 0.0;
    }
}
bool have_IzzMR2() {
    if (IzzMR2.isSet()) {
        return (IzzMR2 > 0.0);
    } else {
        return false;
    }
}

// modified versions of RadioScienceGravityData calls, to include C10,C11,S11 and IzzMR2
unsigned int mod_gravityData_index(const orsaPDS::RadioScienceGravityData * gravityData,
                                   const QString & key) {
    unsigned int index;
    if (key == orsaPDS::RadioScienceGravityData::keyC(0,0)) {
        index = 0;
    } else if (key == orsaPDS::RadioScienceGravityData::keyC(1,0)) {
        index = 1;
    } else if (key == orsaPDS::RadioScienceGravityData::keyC(1,1)) {
        index = 2;
    } else if (key == orsaPDS::RadioScienceGravityData::keyS(1,1)) {
        index = 3;
    } else if (key == "IzzMR2") {
        index = get_IzzMR2_index();
    } else {
        index = gravityData->index(key);
        if (index != 0) index += 3;
    }
    return index;
}

QString mod_gravityData_key(const orsaPDS::RadioScienceGravityData * gravityData,
                            const unsigned int & index) {
    // ORSA_DEBUG("index: %i   gravityData->numberOfCoefficients+3: %i",index,gravityData->numberOfCoefficients);
    if (index == 0) {
        return orsaPDS::RadioScienceGravityData::keyC(0,0);    
    } else if (index == 1) {
        return orsaPDS::RadioScienceGravityData::keyC(1,0);        
    } else if (index == 2) {
        return orsaPDS::RadioScienceGravityData::keyC(1,1);        
    } else if (index == 3) {
        return orsaPDS::RadioScienceGravityData::keyS(1,1);
    } else if (index == get_IzzMR2_index()) {
        return "IzzMR2";
    } else {
        return gravityData->key(index-3);
    }
}
//
double mod_gravityData_getCoeff(const orsaPDS::RadioScienceGravityData * gravityData,
                                const QString & key) {
    double coeff;
    if (key == orsaPDS::RadioScienceGravityData::keyC(0,0)) {
        coeff = 1.0;
    } else if ( (key == orsaPDS::RadioScienceGravityData::keyC(1,0)) ||  
                (key == orsaPDS::RadioScienceGravityData::keyC(1,1)) ||
                (key == orsaPDS::RadioScienceGravityData::keyS(1,1)) ) {
        coeff = 0.0;
    } else if (key == "IzzMR2") {
        coeff = get_IzzMR2();
    } else {
        coeff = gravityData->getCoeff(key);
    }
    return coeff;
}

double mod_gravityData_getCoeff(const orsaPDS::RadioScienceGravityData * gravityData,
                                const unsigned int & index) {
    return mod_gravityData_getCoeff(gravityData,mod_gravityData_key(gravityData,index));
}

unsigned int mod_gravityData_numberOfCoefficients(const orsaPDS::RadioScienceGravityData * gravityData) {
    // return (gravityData->numberOfCoefficients + 3 + (have_IzzMR2() ? 1 : 0));
    unsigned int retVal = gravityData->numberOfCoefficients + 3;
    if (get_IzzMR2_index()==retVal) ++retVal;
    return retVal;
}

gsl_vector * mod_gravityData_getCoefficientVector(const orsaPDS::RadioScienceGravityData * gravityData) {
    gsl_vector * mu = gravityData->getCoefficientVector();
    gsl_vector * mod_mu = gsl_vector_alloc(mod_gravityData_numberOfCoefficients(gravityData));
    for (unsigned int k=0; k<mod_gravityData_numberOfCoefficients(gravityData); ++k) {
        if (k==0) {
            // gsl_vector_set(mod_mu,k,gsl_vector_get(mu,k));
            gsl_vector_set(mod_mu,k,1.0);
        } else if ( (k==1) || (k==2) || (k==3) ) {
            gsl_vector_set(mod_mu,k,0.0);
        } else if (k == get_IzzMR2_index()) {
            gsl_vector_set(mod_mu,k,get_IzzMR2());
        } else {
            gsl_vector_set(mod_mu,k,gsl_vector_get(mu,k-3));
        }
    }
    gsl_vector_free(mu);
    return mod_mu;
}

gsl_matrix * mod_gravityData_getCovarianceMatrix(const orsaPDS::RadioScienceGravityData * gravityData) {
    gsl_matrix * covm = gravityData->getCovarianceMatrix();
    gsl_matrix * mod_covm = gsl_matrix_alloc(mod_gravityData_numberOfCoefficients(gravityData),
                                             mod_gravityData_numberOfCoefficients(gravityData));
    for (unsigned int l=0; l<mod_gravityData_numberOfCoefficients(gravityData); ++l) {
        for (unsigned int m=0; m<mod_gravityData_numberOfCoefficients(gravityData); ++m) { 
            if ((l==0) && (m==0)) {
                // gsl_matrix_set(mod_covm,l,m,gsl_matrix_get(covm,l,m));
                gsl_matrix_set(mod_covm,l,m,gsl_matrix_get(covm,l,m)/orsa::square(gravityData->GM));
#warning NEED TO SCALE THIS DOWN: from GM to 1.0 level... C_{00}... so check this!
            } else if ((l==1) || (l==2) || (l==3) || (l==get_IzzMR2_index()) || (m==1) || (m==2) || (m==3) || (m==get_IzzMR2_index())) {
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
    gsl_matrix_free(covm);
    return mod_covm;   
}

gsl_matrix * mod_gravityData_getInverseCovarianceMatrix(const orsaPDS::RadioScienceGravityData * gravityData) {
    gsl_matrix * inv_covm = gravityData->getInverseCovarianceMatrix();
    gsl_matrix * mod_inv_covm = gsl_matrix_alloc(mod_gravityData_numberOfCoefficients(gravityData),
                                             mod_gravityData_numberOfCoefficients(gravityData));
    for (unsigned int l=0; l<mod_gravityData_numberOfCoefficients(gravityData); ++l) {
        for (unsigned int m=0; m<mod_gravityData_numberOfCoefficients(gravityData); ++m) { 
            if ((l==0) && (m==0)) {
                // gsl_matrix_set(mod_inv_covm,l,m,gsl_matrix_get(inv_covm,l,m));
                gsl_matrix_set(mod_inv_covm,l,m,gsl_matrix_get(inv_covm,l,m)/orsa::square(gravityData->GM));
#warning NEED TO SCALE THIS DOWN: from GM to 1.0 level... C_{00}... so check this!
            } else if ((l==1) || (l==2) || (l==3) || (l==get_IzzMR2_index()) || (m==1) || (m==2) || (m==3) || (m==get_IzzMR2_index())) {
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
        printf("Usage: %s <RadioScienceGravityFile> <plate-model-file> <plate-model-R0_km> <gravity-degree> <polynomial-degree> <CM-x_km> <CM-y_km> <CM-z_km> <CM-sigma-x_km> <CM-sigma-y_km> <CM-sigma-z_km> <Izz/MR^2> <num-sample-points>\n",argv[0]);
        exit(0);
    }   
    
    const double km = orsa::FromUnits(1.0,orsa::Unit::KM);
    
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
    IzzMR2 = atof(argv[12]); // global
    const int numSamplePoints = atoi(argv[13]);
    // const bool have_CCMDF_file = (argc == 15);
    // const std::string CCMDF_filename = (argc == 15) ? argv[14] : "";
    
    // safer over NFS
    // sqlite3_vfs_register(sqlite3_vfs_find("unix-dotfile"), 1);
    
    const std::string SQLiteDBFileName = getSqliteDBFileName_simplex(plateModelFile,plateModelR0);
    
    // test specific cases, for debug purposes only!
    // orsa::GlobalRNG::randomSeed = -800402816;
    
    osg::ref_ptr<orsaPDS::RadioScienceGravityData> gravityData = new orsaPDS::RadioScienceGravityData;
    orsaPDS::RadioScienceGravityFile::read(gravityData.get(),radioScienceGravityFile,512,1518);
    
    osg::ref_ptr<GaskellPlateModel> shapeModel = new GaskellPlateModel;
    if (!shapeModel->read(plateModelFile)) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    osg::ref_ptr<SimplexIntegration<simplex_T> > si = new SimplexIntegration<simplex_T>(shapeModel.get(), plateModelR0, SQLiteDBFileName);
    si->reserve(polynomialDegree);
    
    const GaskellPlateModel::VertexVector & vv = shapeModel->getVertexVector();
    
    const size_t SH_degree = gravityDegree; // shperical harmonics degree
    const size_t  T_degree = polynomialDegree; // chebyshev polynomials degree
    
    // const size_t SH_size = (SH_degree+1)*(SH_degree+1);
    const size_t  SH_size = (SH_degree+1)*(SH_degree+1) + (have_IzzMR2() ? 1 : 0);
    const size_t   T_size = CubicChebyshevMassDistribution::totalSize(T_degree);
    
    IzzMR2_index = (have_IzzMR2() ? SH_size-1 : -1); // global
    
    const double GM = gravityData->GM; 
    
    const double R0_cube_over_M = orsa::cube(plateModelR0)*orsa::Unit::G()/GM;
    
    ORSA_DEBUG("SH_size: %d   T_size: %d   mod_gravityData_numberOfCoefficients: %d",SH_size,T_size,mod_gravityData_numberOfCoefficients(gravityData.get()));
    
    /* 
       if (T_size <= SH_size) {
       ORSA_DEBUG("this method works only when the problem is under-determined, exiting");
       exit(0);
       }
    */
    
    if (gravityDegree > gravityData->degree) {
        ORSA_DEBUG("requested gravity degree [%i] is larger than the RadioScienceGravityFile [%s] degree [%i]; exiting",
                   gravityDegree,radioScienceGravityFile.c_str(),gravityData->degree);
        exit(0);
    }
    
    // osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution;
    /* 
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
    */
    
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

    // SQLite db
    char db_name[4096];
    {
        char arg_of_dirname[4096];
        char arg_of_basename[4096];
        strcpy(arg_of_dirname,SQLiteDBFileName.c_str());
        strcpy(arg_of_basename,SQLiteDBFileName.c_str());
        gmp_sprintf(db_name,"%s/.cT2sh_CM_%.12g_%.12g_%.12g%s",
                    dirname(arg_of_dirname),
                    CMx_over_plateModelR0,
                    CMy_over_plateModelR0,
                    CMz_over_plateModelR0,
                    basename(arg_of_basename));
        // ORSA_DEBUG("db_name: [%s]",db_name);
    }
    sqlite3 * db;
    int rc = sqlite3_open(db_name,&db);
    /* int rc = sqlite3_open_v2(db_name,
       &db,
       SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,
       "unix-dotfile");
    */
    if (rc) {
        fprintf(stderr,"Can't open db: %s\n",sqlite3_errmsg(db));
        sqlite3_close(db);
    }
    
    {
        char * zErr;
        // create results table
        std::string sql = "CREATE TABLE if not exists matrix(z_sh INTEGER, z_cT INTEGER, val REAL, CONSTRAINT un_sh_cT UNIQUE(z_sh,z_cT) ON CONFLICT IGNORE)";
        rc = sqlite3_exec(db,sql.c_str(),NULL,NULL,&zErr);
        //
        if (rc != SQLITE_OK) {
            if (zErr != NULL) {
                fprintf(stderr,"SQL error: %s\n",zErr);
                sqlite3_free(zErr);
            }
        }
    }
    
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
            
            /* const size_t z_C = (l==0) ?
               mod_gravityData_index(gravityData.get(),"GM") :
               mod_gravityData_index(gravityData.get(),orsaPDS::RadioScienceGravityData::keyC(l,m));
            */
            const size_t z_C = mod_gravityData_index(gravityData.get(),orsaPDS::RadioScienceGravityData::keyC(l,m));
            const size_t z_S = (m==0) ? 0 : mod_gravityData_index(gravityData.get(),orsaPDS::RadioScienceGravityData::keyS(l,m));
            
            // ORSA_DEBUG("l=%i m=%i z_C=%i z_S=%i",l,m,z_C,z_S);
            
            {
                // begin transaction
                std::string sql = "begin";
                char * zErr;
                rc = sqlite3_exec(db,sql.c_str(),NULL,NULL,&zErr);
                if (rc != SQLITE_OK) {
                    if (zErr != NULL) {
                        fprintf(stderr,"SQL error: %s\n",zErr);
                        sqlite3_free(zErr);
                    }
                }
            }
            
            // ti,tj,tk are the expansion of the density in terms of the cubic Chebyshev
            for (size_t running_T_degree=0; running_T_degree<=T_degree; ++running_T_degree) {
                for (size_t ti=0; ti<=T_degree; ++ti) {
                    const std::vector<mpz_class> & cTi = orsa::ChebyshevTcoeff(ti);
                    for (size_t tj=0; tj<=T_degree-ti; ++tj) {
                        const std::vector<mpz_class> & cTj = orsa::ChebyshevTcoeff(tj);
                        for (size_t tk=0; tk<=T_degree-ti-tj; ++tk) {
                            
                            if (ti+tj+tk != running_T_degree) continue;
                            
                            const std::vector<mpz_class> & cTk = orsa::ChebyshevTcoeff(tk);
                            
                            const size_t z_cT = CubicChebyshevMassDistribution::index(ti,tj,tk);

                            // ORSA_DEBUG("i=%i j=%i k=%i z_cT=%i",ti,tj,tk,z_cT);
                            
                            // check if available in db
                            
                            bool needToCompute_C=true;
                            {
                                // first check if it is in the SQLite db
                                char **result;
                                int nrows, ncols;
                                char * zErr;
                                char sql_line[1024];
                                sprintf(sql_line,
                                        "SELECT val FROM matrix WHERE z_sh=%zi and z_cT=%zi",
                                        z_C,z_cT);
                                int rc = sqlite3_get_table(db,sql_line,&result,&nrows,&ncols,&zErr);
                                //
                                if (rc != SQLITE_OK) {
                                    if (zErr != NULL) {
                                        fprintf(stderr,"SQL error: %s\n",zErr);
                                        sqlite3_free(zErr);
                                    }
                                }
                                //
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
                                    gsl_matrix_set(cT2sh,z_C,z_cT,atof(result[1])*R0_cube_over_M);
                                    needToCompute_C = false;
                                } else { // if (nrows>size_H) {
                                    ORSA_ERROR("database corrupted, only 1 entry per index is admitted");
                                    ORSA_DEBUG("nrows: %i z_sh=%i z_cT=%i",nrows,z_C,z_cT);
                                    for (int i=0; i<nrows; ++i) {
                                        for (int j=0; j<ncols; ++j) {
                                            // i=0 is the header
                                            const int index = (i+1)*ncols+j;
                                            ORSA_DEBUG("result[%i] = %s",index, result[index]);
                                        }
                                    }
                                }
                                //
                                sqlite3_free_table(result);
                            }
                            
                            bool needToCompute_S=true;
                            if (m==0) needToCompute_S=false;
                            if (m!=0) {
                                // first check if it is in the SQLite db
                                char **result;
                                int nrows, ncols;
                                char * zErr;
                                char sql_line[1024];
                                sprintf(sql_line,
                                        "SELECT val FROM matrix WHERE z_sh=%zi and z_cT=%zi",
                                        z_S,z_cT);
                                int rc = sqlite3_get_table(db,sql_line,&result,&nrows,&ncols,&zErr);
                                //
                                if (rc != SQLITE_OK) {
                                    if (zErr != NULL) {
                                        fprintf(stderr,"SQL error: %s\n",zErr);
                                        sqlite3_free(zErr);
                                    }
                                }
                                // 
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
                                    gsl_matrix_set(cT2sh,z_S,z_cT,atof(result[1])*R0_cube_over_M);
                                    needToCompute_S = false;
                                } else { // if (nrows>size_H) {
                                    ORSA_ERROR("database corrupted, only 1 entry per index is admitted");
                                    ORSA_DEBUG("nrows: %i z_sh=%i z_cT=%i",nrows,z_S,z_cT);
                                    for (int i=0; i<nrows; ++i) {
                                        for (int j=0; j<ncols; ++j) {
                                            // i=0 is the header
                                            const int index = (i+1)*ncols+j;
                                            ORSA_DEBUG("result[%i] = %s",index, result[index]);
                                        }
                                    }
                                }
                                //
                                sqlite3_free_table(result);
                            }
                            
                            if (!needToCompute_C && !needToCompute_S) continue;
                            
                            double C2cT = 0.0;
                            double S2cT = 0.0;
                            
                            // ci,cj,ck are the expansion of each Chebyshev polynomial in terms of powers of x,y,z
                            for (size_t ci=0; ci<=ti; ++ci) {
                                if (cTi[ci] == 0) continue;
                                for (size_t cj=0; cj<=tj; ++cj) {
                                    if (cTj[cj] == 0) continue;
                                    for (size_t ck=0; ck<=tk; ++ck) {
                                        if (cTk[ck] == 0) continue;
                                        // ni,nj,nk are the expansion of C_lm,S_lm in terms of N_ijk
                                        for (size_t ni=0; ni<=l; ++ni) {
                                            for (size_t nj=0; nj<=l-ni; ++nj) {
                                                for (size_t nk=0; nk<=l-ni-nj; ++nk) {
                                                    if ( (C_tri_integral[ni][nj][nk] == 0) && (S_tri_integral[ni][nj][nk] == 0) ) continue;
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
                                    }
                                }
                            }
                            
                            /* C2cT /= si->getIntegral(0,0,0);
                               if (m!=0) S2cT /= si->getIntegral(0,0,0);
                            */
                            
                            // R0^3/M factor
                            C2cT *= R0_cube_over_M;
                            if (m!=0) S2cT *= R0_cube_over_M;
                            
                            // removed the GM factor, now scaling to C_00
                            /* if (l==0) {
                               C2cT *= GM;
                               }
                            */
                            
                            // saving, NOTE how we are adding terms here
                            gsl_matrix_set(cT2sh,z_C,z_cT,gsl_matrix_get(cT2sh,z_C,z_cT)+C2cT);
                            if (m!=0) gsl_matrix_set(cT2sh,z_S,z_cT,gsl_matrix_get(cT2sh,z_S,z_cT)+S2cT);
                            
                            {
                                // save in db
                                char * zErr;
                                char sql_line[1024];
                                sprintf(sql_line,
                                        "INSERT INTO matrix VALUES(%zi,%zi,%.12e)",
                                        z_C,z_cT,gsl_matrix_get(cT2sh,z_C,z_cT)/R0_cube_over_M);
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
                            if (m!=0) {
                                // save in db
                                char * zErr;
                                char sql_line[1024];
                                sprintf(sql_line,
                                        "INSERT INTO matrix VALUES(%zi,%zi,%.12e)",
                                        z_S,z_cT,gsl_matrix_get(cT2sh,z_S,z_cT)/R0_cube_over_M);
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
                }
            }
            
            {
                // close transaction
                std::string sql = "commit";
                char * zErr;
                do {
                    rc = sqlite3_exec(db,sql.c_str(),NULL,NULL,&zErr);
                    if (rc==SQLITE_BUSY) {
                        ORSA_DEBUG("database busy, retrying...");
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
    if (have_IzzMR2()) {
        // one line after all entries above
        // const size_t z_IzzMR2 = 1 + mod_gravityData_index(gravityData.get(),orsaPDS::RadioScienceGravityData::keyS(gravityDegree,gravityDegree));
        const size_t z_IzzMR2 = IzzMR2_index;
        const size_t l = 2;
        
        const double radiusCorrectionFactor = 1.0; // orsa::int_pow(radiusCorrectionRatio,l);
        
        for (size_t running_T_degree=0; running_T_degree<=T_degree; ++running_T_degree) {
            for (size_t ti=0; ti<=T_degree; ++ti) {
                const std::vector<mpz_class> & cTi = orsa::ChebyshevTcoeff(ti);
                for (size_t tj=0; tj<=T_degree-ti; ++tj) {
                    const std::vector<mpz_class> & cTj = orsa::ChebyshevTcoeff(tj);
                    for (size_t tk=0; tk<=T_degree-ti-tj; ++tk) {
                        if (ti+tj+tk != running_T_degree) continue;
                        const std::vector<mpz_class> & cTk = orsa::ChebyshevTcoeff(tk);
                        
                        const size_t z_cT = CubicChebyshevMassDistribution::index(ti,tj,tk);
                        
                        // ORSA_DEBUG("i=%i j=%i k=%i z_cT=%i",ti,tj,tk,z_cT);
                        
                        double IzzMR2_to_cT = 0.0;
                        
                        // ci,cj,ck are the expansion of each Chebyshev polynomial in terms of powers of x,y,z
                        for (size_t ci=0; ci<=ti; ++ci) {
                            if (cTi[ci] == 0) continue;
                            for (size_t cj=0; cj<=tj; ++cj) {
                                if (cTj[cj] == 0) continue;
                                for (size_t ck=0; ck<=tk; ++ck) {
                                    if (cTk[ck] == 0) continue;
                                    // ni,nj,nk are the expansion of IzzMR2 in terms of N_ijk
                                    for (size_t ni=0; ni<=l; ++ni) {
                                        for (size_t nj=0; nj<=l-ni; ++nj) {
                                            for (size_t nk=0; nk<=l-ni-nj; ++nk) {
                                                // if ( (C_tri_integral[ni][nj][nk] == 0) && (S_tri_integral[ni][nj][nk] == 0) ) continue;
                                                // ti,tj,tk are the expansion of the density in terms of the cubic Chebyshev
                                                // IzzMR2
                                                if ( (ni==2 && nj==0 && nk==0) ||
                                                     (ni==0 && nj==2 && nk==0) ) {
                                                    // bi,bj,bk are the binomial expansion about the center of mass
                                                    // this also introduces a power_sign
                                                    for (size_t bi=0; bi<=ni; ++bi) {
                                                        for (size_t bj=0; bj<=nj; ++bj) {
                                                            for (size_t bk=0; bk<=nk; ++bk) {
                                                                
                                                                IzzMR2_to_cT +=
                                                                    orsa::power_sign(bi+bj+bk) *
                                                                    radiusCorrectionFactor *
                                                                    // C_tri_norm[ni][nj][nk] *
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
                            /* C2cT /= si->getIntegral(0,0,0);
                               if (m!=0) S2cT /= si->getIntegral(0,0,0);
                            */
                            
                            // R0^3/M factor
                            IzzMR2_to_cT *= R0_cube_over_M;
                            // if (m!=0) S2cT *= R0_cube_over_M;
                            
                            // removed the GM factor, now scaling to C_00
                            /* if (l==0) {
                               C2cT *= GM;
                               }
                            */
                            
                            // saving, NOTE how we are adding terms here
                            gsl_matrix_set(cT2sh,z_IzzMR2,z_cT,gsl_matrix_get(cT2sh,z_IzzMR2,z_cT)+IzzMR2_to_cT);
                            // if (m!=0) gsl_matrix_set(cT2sh,z_S,z_cT,gsl_matrix_get(cT2sh,z_S,z_cT)+S2cT);
                        }
                    }
                }
            }
        }    
    }
    
    /* for (size_t z_sh=0; z_sh<SH_size; ++z_sh) {
       for (size_t z_cT=0; z_cT<T_size; ++z_cT) {
       size_t Tx,Ty,Tz;
       CubicChebyshevMassDistribution::triIndex(Tx,Ty,Tz,z_cT);
       ORSA_DEBUG("cT2sh[%03i][%03i] = %+9.6f [%s -> cT[%i][%i][%i]]",
       z_sh,z_cT,gsl_matrix_get(cT2sh,z_sh,z_cT),mod_gravityData_key(gravityData.get(),z_sh).toStdString().c_str(),Tx,Ty,Tz);
       }
       }
    */
    
    if (0) {
        // LaTeX output
        ORSA_DEBUG("LaTeX output --------------------");
        // bulk density = M/V = M / (si000*R0^3) 
        const double matrix_factor = 1.0/bulkDensity;
        ORSA_DEBUG(">>>>> NOTE: common factor in front of matrix: %g",matrix_factor);
        for (size_t z_sh=0; z_sh<SH_size; ++z_sh) {
            for (size_t z_cT=0; z_cT<T_size; ++z_cT) {
                const double val = gsl_matrix_get(cT2sh,z_sh,z_cT) / matrix_factor;
                gmp_fprintf(stderr,"%9.6f & ",val);
                // gmp_fprintf(stderr,"@%02i/%02i %9.6f &",z_sh,z_cT,val);
            }   
            gmp_fprintf(stderr,"\n");
        }
        ORSA_DEBUG("---------------------------------");
    }
    
    gsl_vector * pds_coeff    = mod_gravityData_getCoefficientVector(gravityData.get());
    gsl_matrix * pds_covm     = mod_gravityData_getCovarianceMatrix(gravityData.get());
    // gsl_matrix * pds_inv_covm = mod_gravityData_getInverseCovarianceMatrix(gravityData.get());
    
    gsl_matrix_view     pds_covm_view = gsl_matrix_submatrix(    pds_covm, 0, 0, SH_size, SH_size);
    // gsl_matrix_view pds_inv_covm_view = gsl_matrix_submatrix(pds_inv_covm, 0, 0, SH_size, SH_size);
    
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
        
        ORSA_DEBUG("null space size: %i",N-M);
        
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
                    ORSA_DEBUG("uK[%03i][%03i] =%12.6f (null space base vector for cT[%i][%i][%i])",b,s,gsl_vector_get(uK[b],s),Tx,Ty,Tz);
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
        
        if (0) {
            // LaTeX output
            ORSA_DEBUG("LaTeX output --------------------");
            const double matrix_factor = bulkDensity;
            ORSA_DEBUG(">>>>> NOTE: common factor in front of matrix: %g",matrix_factor);
            // const double eps = 1.0e-12; // values smaller than eps are printed as zeros
            for (size_t row=0; row<N; ++row) {
                for (size_t col=0; col<M; ++col) {
                    const double val = gsl_matrix_get(pseudoInvA,row,col) / matrix_factor;
                    gmp_fprintf(stderr,"%12.6f & ",val);
                }
                gmp_fprintf(stderr,"\n");
            }
            ORSA_DEBUG("---------------------------------");
        }
        
        
        gsl_vector * sh = gsl_vector_calloc(M); // b  of A x = b
        gsl_vector * cT = gsl_vector_calloc(N); // x  of A x = b

        // gsl_vector * alt_sh = gsl_vector_calloc(M);
        
        // used in the loop...
        gsl_vector * cT0 = gsl_vector_calloc(N);
        
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
                ORSA_ERROR("problems with the covariance matrix: null eigenvalue found, i=%i",i);
            }
            sigma[i] = sqrt(fabs(gsl_vector_get(eval,i)));
            // ORSA_DEBUG("sigma[%i] = %g",i,sigma[i]);
        }
        gsl_vector * sampleCoeff_x  = gsl_vector_alloc(M);
        gsl_vector * sampleCoeff_y  = gsl_vector_alloc(M); 
        
        // const double minimumDensity = orsa::FromUnits(orsa::FromUnits(1.50,orsa::Unit::GRAM),orsa::Unit::CM,-3);
        
        // output file, global
        {
            char line[4096];
            gmp_sprintf(line,
                        // "CCMDF_%i_%i_%g_%i.out",
                        "CCMDF_%i_%i_%i.out",
                        gravityDegree,
                        polynomialDegree,
                        // orsa::FromUnits(orsa::FromUnits(minimumDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
                        getpid());
            CCMDF_output_filename = line;
        }
        
        // for (size_t gen=0; gen<1000; ++gen) {
        // for (size_t gen=0; gen<10; ++gen) {
        while (1) {
            
            // osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution;
            
            // set massDistribution
            CubicChebyshevMassDistribution::CoefficientType coeff;
            CubicChebyshevMassDistribution::resize(coeff,T_degree);
            for (size_t runningDegree=0; runningDegree<=T_degree; ++runningDegree) {
                for (size_t i=0; i<=T_degree; ++i) {
                    for (size_t j=0; j<=T_degree-i; ++j) {
                        for (size_t k=0; k<=T_degree-i-j; ++k) {
                            if (i+j+k == runningDegree) {
                                coeff[i][j][k] = 0.0;
                                // ORSA_DEBUG("ijk: %i %i %i   coeff: %g",i,j,k,data.coeff[i][j][k]);
                            }
                        }
                    }
                }
            }
            // set coeff above, or just set the non-zero coeff here
            coeff[0][0][0] = bulkDensity;
            
            /* for (size_t runningDegree=0; runningDegree<=T_degree; ++runningDegree) {
               for (size_t i=0; i<=T_degree; ++i) {
               for (size_t j=0; j<=T_degree-i; ++j) {
               for (size_t k=0; k<=T_degree-i-j; ++k) {
               if (i+j+k == runningDegree) {
               // coeff[i][j][k] = 0.0;
               ORSA_DEBUG("ijk: %i %i %i   coeff: %g",i,j,k,coeff[i][j][k]);
               }
               }
               }
               }
               }
            */
            
            std::vector< std::vector<mpf_class> > NOLayerData_norm_C;
            std::vector< std::vector<mpf_class> > NOLayerData_norm_S;
            //
            NOLayerData_norm_C.resize(SH_degree+1);
            NOLayerData_norm_S.resize(SH_degree+1);
            for (size_t l=0; l<=SH_degree; ++l) {
                NOLayerData_norm_C[l].resize(l+1);
                NOLayerData_norm_S[l].resize(l+1);
                for (size_t m=0; m<=l; ++m) {
                    NOLayerData_norm_C[l][m] = 0.0;
                    NOLayerData_norm_S[l][m] = 0.0;
                }
            }
            mpf_class NOLayerData_IzzMR2 = 0.0;
            {
                osg::ref_ptr<CubicChebyshevMassDistribution> md_NOLD =
                    new CubicChebyshevMassDistribution(coeff,
                                                       plateModelR0,
                                                       0);
                orsa::Cache<orsa::Vector> CM = sampled_CM;
                CM.lock();
                CCMD2SH(CM,
                        NOLayerData_norm_C,
                        NOLayerData_norm_S,
                        NOLayerData_IzzMR2,
                        SH_degree, // gravityData->degree,
                        si.get(),
                        md_NOLD,
                        plateModelR0,
                        gravityData->R0);
                
                // scale later, after massDistribution has been created
                /* const double NOLayerMassFraction = 1.0 - massDistribution->layerData->totalExcessMass() / (GM/orsa::Unit::G());
                   for (size_t l=0; l<=SH_degree; ++l) {
                   for (size_t m=0; m<=l; ++m) {
                   NOLayerData_norm_C[l][m] *= NOLayerMassFraction;
                   NOLayerData_norm_S[l][m] *= NOLayerMassFraction;
                   }
                   }
                   NOLayerData_IzzMR2 *= NOLayerMassFraction;
                */
            }

            if (1) {
                // test
                std::vector<double> gravity_data;
                std::vector<double> gravity_CCMD;
                gravity_data.reserve(2*SH_degree+1);
                gravity_CCMD.reserve(2*SH_degree+1);
                for (size_t l=0; l<=SH_degree; ++l) {
                    gravity_data.clear();
                    gravity_CCMD.clear();
                    for (size_t m=1; m<=l; ++m) {
                        const QString key = orsaPDS::RadioScienceGravityData::keyS(l,m);
                        gravity_data.push_back(mod_gravityData_getCoeff(gravityData.get(),key));
                        gravity_CCMD.push_back(NOLayerData_norm_S[l][m].get_d());
                    }
                    for (size_t m=0; m<=l; ++m) {
                        const QString key = orsaPDS::RadioScienceGravityData::keyC(l,m);
                        gravity_data.push_back(mod_gravityData_getCoeff(gravityData.get(),key));
                        gravity_CCMD.push_back(NOLayerData_norm_C[l][m].get_d());
                    }
                    double scalar_product=0.0;
                    double l2_CCMD=0.0;
                    for (size_t k=0; k<gravity_data.size(); ++k) {
                        scalar_product += gravity_data[k]*gravity_CCMD[k];
                        l2_CCMD += orsa::square(gravity_CCMD[k]);
                        // ORSA_DEBUG("%g %g",gravity_data[k],gravity_CCMD[k]);
                    }
                    const double factor = scalar_product/l2_CCMD;
                    ORSA_DEBUG("l = %2i  factor: %g",l,factor);
                    
                }
            }            
            
            // FastLayersMultifit 
            osg::ref_ptr<orsa::MultifitParameters> par = new orsa::MultifitParameters;
            std::vector<FastLayersMultifit::EllipsoidLayerData> ellipsoidLayerData;
            {
                ellipsoidLayerData.resize(2); // number of ellipsoid layers

                {
                    double totalExcessMass;
                    do {
                        
                        for (size_t k=0; k<ellipsoidLayerData.size(); ++k) {
                            ellipsoidLayerData[k].excessMass = (GM/orsa::Unit::G())*(0.00+1.00*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform());
                            
                        }
                        
                        totalExcessMass=0.0;
                        for (size_t k=0; k<ellipsoidLayerData.size(); ++k) {
                            totalExcessMass += ellipsoidLayerData[k].excessMass;
                        }
                        
                    } while (totalExcessMass >= (GM/orsa::Unit::G()));
                }
                
                char varName[4096];
                
                // same for all
                sprintf(varName,"v0x");
                par->insert(varName,(0.010+0.020*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform())*km,0.001*km);
                sprintf(varName,"v0y");
                par->insert(varName,(0.010+0.020*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform())*km,0.001*km);
                sprintf(varName,"v0z");
                par->insert(varName,(0.010+0.020*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform())*km,0.001*km);
                //
                sprintf(varName,"psi");
                par->insert(varName,orsa::twopi()*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform(),0.00001);
                sprintf(varName,"theta");
                par->insert(varName,orsa::twopi()*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform(),0.00001);
                sprintf(varName,"phi");
                par->insert(varName,orsa::twopi()*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform(),0.00001);
                
                for (size_t k=0; k<ellipsoidLayerData.size(); ++k) {
                    
                    ORSA_DEBUG("ellipsoidLayerData[%i] mass fraction: %g",k,ellipsoidLayerData[k].excessMass/(GM/orsa::Unit::G()));
                    
                    sprintf(varName,"a_%i",k);
                    par->insert(varName,(80.0+40.0*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform())*km,0.001*km);
                    // par->setRange(varName,0.0*km,300.0*km);
                    sprintf(varName,"b_%i",k);
                    par->insert(varName,(80.0+40.0*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform())*km,0.001*km);
                    // par->setRange(varName,0.0*km,300.0*km);
                    sprintf(varName,"c_%i",k);
                    par->insert(varName,(80.0+40.0*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform())*km,0.001*km);
                    // par->setRange(varName,0.0*km,300.0*km);
                    //
                    /* sprintf(varName,"v0x_%i",k);
                       par->insert(varName,(0.010+0.020*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform())*km,0.001*km);
                       sprintf(varName,"v0y_%i",k);
                       par->insert(varName,(0.010+0.020*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform())*km,0.001*km);
                       sprintf(varName,"v0z_%i",k);
                       par->insert(varName,(0.010+0.020*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform())*km,0.001*km);
                       //
                       sprintf(varName,"psi_%i",k);
                       par->insert(varName,orsa::twopi()*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform(),0.00001);
                       sprintf(varName,"theta_%i",k);
                       par->insert(varName,orsa::twopi()*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform(),0.00001);
                       sprintf(varName,"phi_%i",k);
                       par->insert(varName,orsa::twopi()*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform(),0.00001);
                    */
                }
            }
            
            osg::ref_ptr<FastLayersMultifit> mf = new FastLayersMultifit;
            mf->ellipsoidLayerData = ellipsoidLayerData;
            mf->totalBodyMass = (GM/orsa::Unit::G());
            mf->SH_degree = 4;
            mf->SH_size = (mf->SH_degree+1)*(mf->SH_degree+1);
            mf->R0_plate = plateModelR0;
            mf->R0_gravity = gravityData->R0;
            mf->sampled_CM = sampled_CM;
            mf->si = si;
            mf->gravityData = gravityData;
            
            {
                // scale NOLayerData SH
                double ellipsoidTotalExcessMass = 0.0;
                for (size_t k=0; k<ellipsoidLayerData.size(); ++k) {
                    ellipsoidTotalExcessMass += ellipsoidLayerData[k].excessMass;
                }
                const double factor = 1.0 - ellipsoidTotalExcessMass/(GM/orsa::Unit::G());

                for (size_t l=0; l<=SH_degree; ++l) {
                    for (size_t m=0; m<=l; ++m) {
                        NOLayerData_norm_C[l][m] *= factor;
                        if (m!=0) NOLayerData_norm_S[l][m] *= factor;
                    }
                }
                
                /* for (size_t l=0; l<=SH_degree; ++l) {
                   for (size_t m=0; m<=l; ++m) {
                   ORSA_DEBUG("NOLayerData_norm_C[%i][%i] = %g",l,m,NOLayerData_norm_C[l][m].get_d());
                   if (m!=0) ORSA_DEBUG("NOLayerData_norm_S[%i][%i] = %g",l,m,NOLayerData_norm_S[l][m].get_d());
                   }
                   }
                */
            }
            
            osg::ref_ptr<orsa::MultifitData> data = new orsa::MultifitData;
            {
                gsl_matrix * pds_covm  = mod_gravityData_getCovarianceMatrix(gravityData.get());
                for (size_t row=0; row<mf->SH_size; ++row) {
                    orsa::Cache<double> NO_layer_coeff;
                    const QString ref_key = mod_gravityData_key(gravityData.get(),row);
                    for (size_t l=0; l<=SH_degree; ++l) {
                        for (size_t m=0; m<=l; ++m) {
                            if (orsaPDS::RadioScienceGravityData::keyC(l,m) == ref_key) {
                                // ORSA_DEBUG("found: [%s] for l=%i, m=%i",ref_key.toStdString().c_str(),l,m);
                                // layer_coeff = layerData_norm_C[l][m].get_d();
                                NO_layer_coeff = NOLayerData_norm_C[l][m].get_d();
                            } else if (orsaPDS::RadioScienceGravityData::keyS(l,m) == ref_key) {
                                // ORSA_DEBUG("found: [%s] for l=%i, m=%i",ref_key.toStdString().c_str(),l,m);
                                // layer_coeff = layerData_norm_S[l][m].get_d();
                                NO_layer_coeff = NOLayerData_norm_S[l][m].get_d();
                            }                       
                        }
                    }
                    if (!NO_layer_coeff.isSet()) {
                        ORSA_DEBUG("problems... M=%i",M);
                        exit(0);
                    }
                    // ORSA_DEBUG("NO_layer_coeff = %g",(*NO_layer_coeff));
                    // data->insertF(row,gsl_vector_get(alt_sh,row));
                    data->insertF(row,gsl_vector_get(pds_coeff,row)-NO_layer_coeff);
                    // double sigma = sqrt(gsl_matrix_get(pds_covm,row,row));
                    // if (sigma == 0.0) sigma = 1.0e-3*data->getF(row); // for degree 1 terms, which have zero covariance
                    // this sigma is not the one from data, but it's arbitrarily chosen to make the algorithm converge
                    // const double sigma = 1.0e-6 + 0.01*fabs(data->getF(row));
                    const double sigma = 1.0e-9;
                    data->insertSigma(row,sigma);
                    /* ORSA_DEBUG("F for index %03i = [%7s] = %12.6g +/- %12.6g",
                       row,
                       mod_gravityData_key(gravityData.get(),row).toStdString().c_str(),
                       data->getF(row),
                       data->getSigma(row));
                    */
                }
                gsl_matrix_free(pds_covm);
            }
            
            mf->setMultifitParameters(par.get());
            mf->setMultifitData(data.get());
            // mf->setLogFile("FastLayersMultifit.log");
            const bool success = mf->run();
            if (!success) continue;
            
            /* for (size_t row=0; row<M; ++row) {
               const double f = data->getF(row);
               const double T = mf->fun(par.get(),
               data.get(),
               0,
               0,
               row);
               const double err = T-f;
               ORSA_DEBUG("FINAL: %12.6g %12.6g %12.6g",
               f,T,err);
               }
            */
            
            for (size_t s=0; s<par->totalSize(); ++s) {
                ORSA_DEBUG("par[%03i] = [%s] = %12.6g",s,par->name(s).c_str(),par->get(s));
            }
            
            LayerData::EllipsoidLayerVectorType ellipsoidLayerVector;
            {
                bool skip=false;
                char varName[4096];
                
                sprintf(varName,"v0x");
                const double v0x = par->get(varName);
                sprintf(varName,"v0y");
                const double v0y = par->get(varName);
                sprintf(varName,"v0z");
                const double v0z = par->get(varName);
                //
                sprintf(varName,"psi");
                const double psi = par->get(varName);
                sprintf(varName,"theta");
                const double theta = par->get(varName);
                sprintf(varName,"phi");
                const double phi = par->get(varName);
                
                for (size_t k=0; k<ellipsoidLayerData.size(); ++k) {
                    sprintf(varName,"a_%i",k);
                    const double a = fabs(par->get(varName));
                    sprintf(varName,"b_%i",k);
                    const double b = fabs(par->get(varName));
                    sprintf(varName,"c_%i",k);
                    const double c = fabs(par->get(varName));
                    //
                    /* sprintf(varName,"v0x_%i",k);
                       const double v0x = par->get(varName);
                       sprintf(varName,"v0y_%i",k);
                       const double v0y = par->get(varName);
                       sprintf(varName,"v0z_%i",k);
                       const double v0z = par->get(varName);
                       //
                       sprintf(varName,"psi_%i",k);
                       const double psi = par->get(varName);
                       sprintf(varName,"theta_%i",k);
                       const double theta = par->get(varName);
                       sprintf(varName,"phi_%i",k);
                       const double phi = par->get(varName);
                    */
                    
                    const double volume = 4.0/3.0*orsa::pi()*a*b*c;
                    const double excessDensity = ellipsoidLayerData[k].excessMass/volume;
                    
                    ORSA_DEBUG("excessDensity[%i] = %g",k,excessDensity);
                    
                    orsa::Matrix rot;
                    orsa::eulerAnglesToMatrix(rot,
                                              psi,
                                              theta,
                                              phi);
                    
                    osg::ref_ptr<LayerData::EllipsoidLayer> ellipsoidLayer = new LayerData::EllipsoidLayer(excessDensity,a,b,c,orsa::Vector(v0x,v0y,v0z),rot);
                    
                    // first, quicker test
                    if (std::max(a,std::max(b,c))>shapeModel->boundingRadius()) {
                        ORSA_DEBUG("===> layer not fully contained in shape, skipping... *****");
                        skip=true;
                        break;
                    }
                    if (skip) break;
                    // detailed test...
                    for (size_t p=0; p<vv.size(); ++p) {
                        if (ellipsoidLayer->containsPoint(vv[p])) {
                            ORSA_DEBUG("===> layer not fully contained in shape, skipping... *****");
                            skip=true;
                            break;
                        }
                    }
                    if (skip) break;
                    
                    // ellipsoidLayerVector.push_back(new LayerData::EllipsoidLayer(excessDensity,a,b,c,orsa::Vector(v0x,v0y,v0z),rot));
                    ellipsoidLayerVector.push_back(ellipsoidLayer);
                }
                if (skip) continue;
            }
            LayerData::SHLayerVectorType shLayerVector;
            osg::ref_ptr<const LayerData> layerData = new LayerData(ellipsoidLayerVector,shLayerVector);
            {
                // before creating the massDistribution, scale the coeff in order to conserve total mass
                const double coeff_target_mass = (GM/orsa::Unit::G()) - layerData->totalExcessMass();
                std::vector< std::vector< std::vector<double> > > N;
                osg::ref_ptr<CubicChebyshevMassDistribution> tmp_MD = new CubicChebyshevMassDistribution(coeff,plateModelR0,0);
                CCMD2ijk(N,
                         0,
                         si.get(),
                         tmp_MD.get(),
                         plateModelR0);
                const double coeff_current_mass = N[0][0][0]*orsa::cube(plateModelR0);
                const double coeff_factor = coeff_target_mass/coeff_current_mass;
                for (size_t i=0; i<=T_degree; ++i) {
                    for (size_t j=0; j<=T_degree-i; ++j) {
                        for (size_t k=0; k<=T_degree-i-j; ++k) {
                            coeff[i][j][k] *= coeff_factor;
                        }
                    }
                }
            }
            osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
                new CubicChebyshevMassDistribution(coeff,
                                                   plateModelR0,
                                                   layerData);
            
            std::vector< std::vector<mpf_class> > layerData_norm_C;
            std::vector< std::vector<mpf_class> > layerData_norm_S;
            //
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
            mpf_class layerData_IzzMR2 = 0.0;
            if (massDistribution.get() != 0) {
                if (massDistribution->layerData.get() != 0) {
                    CubicChebyshevMassDistribution::CoefficientType md_lD_coeff;
                    CubicChebyshevMassDistribution::resize(md_lD_coeff,0);
                    md_lD_coeff[0][0][0] = 0;
                    osg::ref_ptr<CubicChebyshevMassDistribution> md_lD =
                        new CubicChebyshevMassDistribution(md_lD_coeff,
                                                           plateModelR0,
                                                           massDistribution->layerData.get());
                    orsa::Cache<orsa::Vector> CM = sampled_CM;
                    CM.lock();
                    CCMD2SH(CM,
                            layerData_norm_C,
                            layerData_norm_S,
                            layerData_IzzMR2,
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
                    layerData_IzzMR2 *= layerMassFraction;
                }
            }
            
            // std::vector< std::vector<mpf_class> > NOLayerData_norm_C;
            // std::vector< std::vector<mpf_class> > NOLayerData_norm_S;
            //
            NOLayerData_norm_C.resize(SH_degree+1);
            NOLayerData_norm_S.resize(SH_degree+1);
            for (size_t l=0; l<=SH_degree; ++l) {
                NOLayerData_norm_C[l].resize(l+1);
                NOLayerData_norm_S[l].resize(l+1);
                for (size_t m=0; m<=l; ++m) {
                    NOLayerData_norm_C[l][m] = 0.0;
                    NOLayerData_norm_S[l][m] = 0.0;
                }
            }
            NOLayerData_IzzMR2 = 0.0;
            if (massDistribution.get() != 0) {
                osg::ref_ptr<CubicChebyshevMassDistribution> md_NOLD =
                    new CubicChebyshevMassDistribution(massDistribution->coeff, 
                                                       plateModelR0,
                                                       0);
                orsa::Cache<orsa::Vector> CM = sampled_CM;
                CM.lock();
                CCMD2SH(CM,
                        NOLayerData_norm_C,
                        NOLayerData_norm_S,
                        NOLayerData_IzzMR2,
                        SH_degree, // gravityData->degree,
                        si.get(),
                        md_NOLD,
                        plateModelR0,
                        gravityData->R0);
                
                const double NOLayerMassFraction = 1.0 - massDistribution->layerData->totalExcessMass() / (GM/orsa::Unit::G());
                for (size_t l=0; l<=SH_degree; ++l) {
                    for (size_t m=0; m<=l; ++m) {
                        NOLayerData_norm_C[l][m] *= NOLayerMassFraction;
                        NOLayerData_norm_S[l][m] *= NOLayerMassFraction;
                    }
                }
                NOLayerData_IzzMR2 *= NOLayerMassFraction;
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
                orsa::Cache<double> NO_layer_coeff;
                if (massDistribution.get() == 0) {
                    layer_coeff = 0.0;
                    NO_layer_coeff = 0.0;
                } else if (massDistribution->layerData.get() == 0) {
                    layer_coeff = 0.0;
                } else {
                    const QString ref_key = mod_gravityData_key(gravityData.get(),i);
                    // if (ref_key == "GM") {
                    /* if (ref_key == orsaPDS::RadioScienceGravityData::keyC(0,0)) {
                    // ORSA_DEBUG("found: [%s]",ref_key.toStdString().c_str());
                    // layer_coeff = massDistribution->layerData->totalExcessMass()*orsa::Unit::G();
                    } else {
                    */
                    for (size_t l=0; l<=SH_degree; ++l) {
                        for (size_t m=0; m<=l; ++m) {
                            if (orsaPDS::RadioScienceGravityData::keyC(l,m) == ref_key) {
                                // ORSA_DEBUG("found: [%s] for l=%i, m=%i",ref_key.toStdString().c_str(),l,m);
                                layer_coeff = layerData_norm_C[l][m].get_d();
                                NO_layer_coeff = NOLayerData_norm_C[l][m].get_d();
                            } else if (orsaPDS::RadioScienceGravityData::keyS(l,m) == ref_key) {
                                // ORSA_DEBUG("found: [%s] for l=%i, m=%i",ref_key.toStdString().c_str(),l,m);
                                layer_coeff = layerData_norm_S[l][m].get_d();
                                NO_layer_coeff = NOLayerData_norm_S[l][m].get_d();
                            }                       
                        }
                    }
                    if ("IzzMR2" == ref_key) {
                        layer_coeff = layerData_IzzMR2.get_d();
                        NO_layer_coeff = NOLayerData_IzzMR2.get_d();
                    }
                }
                if (!layer_coeff.isSet()) {
                    ORSA_DEBUG("problems... M=%i",M);
                    exit(0);
                }
                if (!NO_layer_coeff.isSet()) {
                    ORSA_DEBUG("problems... M=%i",M);
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
                
                /* ORSA_DEBUG("%7s = %12.6g [sampled] = %12.6g [layers] + %12.6g   nominal: %+12.6g   delta: %+12.6g   sigma: %12.6g",
                   mod_gravityData_key(gravityData.get(),i).toStdString().c_str(),
                   gsl_vector_get(sh,i)+(*layer_coeff),
                   (*layer_coeff),
                   gsl_vector_get(sh,i),
                   mod_gravityData_getCoeff(gravityData.get(),mod_gravityData_key(gravityData.get(),i)),
                   gsl_vector_get(sh,i)+(*layer_coeff)-mod_gravityData_getCoeff(gravityData.get(),mod_gravityData_key(gravityData.get(),i)),
                   sqrt(gsl_matrix_get(pds_covm,i,i)));
                */
                
                ORSA_DEBUG("%7s = %12.6g [sampled] = %12.6g [CCMD] + %12.6g [layers] + %12.6g   nominal: %+12.6g   delta: %+12.6g   sigma: %12.6g",
                           mod_gravityData_key(gravityData.get(),i).toStdString().c_str(),
                           gsl_vector_get(sh,i)+(*layer_coeff),
                           (*NO_layer_coeff),
                           (*layer_coeff),
                           gsl_vector_get(sh,i)-(*NO_layer_coeff),
                           mod_gravityData_getCoeff(gravityData.get(),mod_gravityData_key(gravityData.get(),i)),
                           gsl_vector_get(sh,i)+(*layer_coeff)-mod_gravityData_getCoeff(gravityData.get(),mod_gravityData_key(gravityData.get(),i)),
                           sqrt(gsl_matrix_get(pds_covm,i,i)));
                
                gsl_matrix_free(pds_covm);
            }
            
            // solving here!
            gsl_blas_dgemv(CblasNoTrans,1.0,pseudoInvA,sh,0.0,cT);

            // allocated earlier...
            // gsl_vector * cT0 = gsl_vector_calloc(N);
            gsl_vector_memcpy(cT0,cT);

            if (0) {
                for (size_t s=0; s<N; ++s) {
                    ORSA_DEBUG("cT0[%i] = %g",s,gsl_vector_get(cT0,s));
                }
            }

            if (0) {
                for (size_t s=0; s<N; ++s) {
                    ORSA_DEBUG("cT0[%i] = %12.6f * bulkDensity",s,gsl_vector_get(cT0,s)/bulkDensity);
                }
            }
            
            if (1) {
                
                // just project on the solution space
                
                const size_t uK_size = N-M;
                std::vector<double> factor;
                factor.resize(uK_size);
                {
                    const size_t cT_CCMD_degree = massDistribution->coeff.size()-1;
                    const size_t cT_CCMD_size   = CubicChebyshevMassDistribution::totalSize(cT_CCMD_degree);
                    size_t Tx,Ty,Tz;
                    for (size_t b=0; b<uK_size; ++b) {
                        factor[b] = 0.0;
                        for (size_t s=0; s<N; ++s) {
                            double cT_CCMD = 0.0; // default value
                            if (s < cT_CCMD_size) {
                                CubicChebyshevMassDistribution::triIndex(Tx,Ty,Tz,s);
                                cT_CCMD  = massDistribution->coeff[Tx][Ty][Tz];
                            }
                            factor[b] += (cT_CCMD-gsl_vector_get(cT0,s))*gsl_vector_get(uK[b],s);
                        }
                        // ORSA_DEBUG("factor[%03i] = %g",b,factor[b]);
                    }
                }
                
                gsl_vector_memcpy(cT,cT0);
                for (size_t b=0; b<uK_size; ++b) {
                    for (size_t j=0; j<T_size; ++j) {
                        gsl_vector_set(cT,j,gsl_vector_get(cT0,j)+factor[b]*gsl_vector_get(uK[b],j));
                    }
                }
                
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
                
                osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
                    new CubicChebyshevMassDistribution(coeff,
                                                       plateModelR0,
                                                       layerData);
                
                {
                    // another quick output...
                    CubicChebyshevMassDistributionFile::CCMDF_data data;
                    data.minDensity = 0.0;
                    data.maxDensity = 0.0;
                    data.deltaDensity = 0.0;
                    data.penalty = 0.0;
                    // data.densityScale = x->bulkDensity;
                    data.R0 = plateModelR0;
                    data.SH_degree = SH_degree;
                    data.coeff = coeff;
                    data.layerData = layerData;
                    //
                    if (1) {
                        orsa::Cache<orsa::Vector> CM = sampled_CM;
                        mpf_class IxxMR2, IyyMR2, IzzMR2;
                        inertia(CM,
                                IxxMR2,
                                IyyMR2,
                                IzzMR2,
                                // SH_degree,
                                si.get(),
                                massDistribution.get(),
                                plateModelR0);
                        
                        char comment[4096];
                        sprintf(comment,"%.6f %.6f %.6f ",IxxMR2.get_d(),IyyMR2.get_d(),IzzMR2.get_d());
                        if (0) {
                            char tmpstr[4096];
                            for (size_t b=0; b<uK_size; ++b) {
                                sprintf(tmpstr," %+12.6f",factor[b]);
                                strcat(comment,tmpstr);
                            }
                        }
                        if (1) {
                            char tmpstr[4096];
                            for (size_t k=0; k<ellipsoidLayerData.size(); ++k) {
                                sprintf(tmpstr," %+12.6f",ellipsoidLayerData[k].excessMass/(GM/orsa::Unit::G()));
                                strcat(comment,tmpstr);
                            }
                        }
                        // strcat(comment," ");
                        // strcat(comment,pvline);
                        data.comment = comment;
                    }
                    //
                    // CubicChebyshevMassDistributionFile::append(data,"CCMDF.out");
                    CubicChebyshevMassDistributionFile::append(data,CCMDF_output_filename.c_str());
                }
                
                
                
            }
            
        }
        
        // free GSL stuff
        gsl_matrix_free(A);
        gsl_matrix_free(AT);
        gsl_matrix_free(QR);
        gsl_vector_free(tau);
        gsl_matrix_free(Q);
        gsl_matrix_free(R);
        for (size_t b=0; b<(N-M); ++b) {
            gsl_vector_free(uK[b]);
        }
        gsl_matrix_free(A_AT);
        gsl_matrix_free(inv_A_AT);
        gsl_matrix_free(pseudoInvA);
        gsl_vector_free(sh);
        // gsl_vector_free(alt_sh);
        gsl_vector_free(cT);
        gsl_eigen_symmv_free(w);
        gsl_vector_free(eval);
        gsl_matrix_free(evec);
        gsl_vector_free(sampleCoeff_x);
        gsl_vector_free(sampleCoeff_y);
        gsl_vector_free(cT0);
        
    }
    
    // free GSL stuff
    gsl_vector_free(pds_coeff);
    gsl_matrix_free(pds_covm);
    // gsl_matrix_free(pds_inv_covm);
    gsl_matrix_free(cT2sh);

    sqlite3_close(db);
    
    return 0;
}

// multifit

double FastLayersMultifit::fun(const orsa::MultifitParameters * par, 
                               const orsa::MultifitData       * data,
                               const unsigned int p, 
                               const int          d,
                               const unsigned int row) const {
    char varName[4096];
    
    LayerData::EllipsoidLayerVectorType ellipsoidLayerVector;
    
    sprintf(varName,"v0x");
    double v0x = par->get(varName);
    if (p == par->index(varName))
        v0x += d*par->getDelta(varName);
    sprintf(varName,"v0y");
    double v0y = par->get(varName);
    if (p == par->index(varName))
        v0y += d*par->getDelta(varName);
    sprintf(varName,"v0z");
    double v0z = par->get(varName);
    if (p == par->index(varName))
        v0z += d*par->getDelta(varName);
    
    sprintf(varName,"psi");
    double psi   = par->get(varName);
    if (p == par->index(varName))
        psi += d*par->getDelta(varName);
    sprintf(varName,"theta");
    double theta = par->get(varName);
    if (p == par->index(varName))
        theta += d*par->getDelta(varName);
    sprintf(varName,"phi");
    double phi   = par->get(varName);
    if (p == par->index(varName))
        phi += d*par->getDelta(varName);
    
    for (size_t k=0; k<ellipsoidLayerData.size(); ++k) {
        
        sprintf(varName,"a_%i",k);
        double a = par->get(varName);
        if (p == par->index(varName))
            a += d*par->getDelta(varName);
        sprintf(varName,"b_%i",k);
        double b = par->get(varName);
        if (p == par->index(varName))
            b += d*par->getDelta(varName);
        sprintf(varName,"c_%i",k);
        double c = par->get(varName);
        if (p == par->index(varName))
            c += d*par->getDelta(varName);
        
        a = fabs(a);
        b = fabs(b);
        c = fabs(c);
        
        /* sprintf(varName,"v0x_%i",k);
           double v0x = par->get(varName);
           if (p == par->index(varName))
           v0x += d*par->getDelta(varName);
           sprintf(varName,"v0y_%i",k);
           double v0y = par->get(varName);
           if (p == par->index(varName))
           v0y += d*par->getDelta(varName);
           sprintf(varName,"v0z_%i",k);
           double v0z = par->get(varName);
           if (p == par->index(varName))
           v0z += d*par->getDelta(varName);
           
           sprintf(varName,"psi_%i",k);
           double psi   = par->get(varName);
           if (p == par->index(varName))
           psi += d*par->getDelta(varName);
           sprintf(varName,"theta_%i",k);
           double theta = par->get(varName);
           if (p == par->index(varName))
           theta += d*par->getDelta(varName);
           sprintf(varName,"phi_%i",k);
           double phi   = par->get(varName);
           if (p == par->index(varName))
           phi += d*par->getDelta(varName);
        */
        
        orsa::Matrix rot;
        eulerAnglesToMatrix(rot,
                            psi,
                            theta,
                            phi);
        
        const double volume = 4.0/3.0*orsa::pi()*a*b*c;
        const double excessDensity = ellipsoidLayerData[k].excessMass/volume;
        
        ellipsoidLayerVector.push_back(new LayerData::EllipsoidLayer(excessDensity,
                                                                     a,b,c,
                                                                     orsa::Vector(v0x,v0y,v0z),
                                                                     rot));
    }
    
    LayerData::SHLayerVectorType emptySHLayerVector;
    osg::ref_ptr<LayerData> layerData = new LayerData(ellipsoidLayerVector,
                                                      emptySHLayerVector);
    
    std::vector< std::vector<mpf_class> > layerData_norm_C;
    std::vector< std::vector<mpf_class> > layerData_norm_S;
    mpf_class layerData_IzzMR2;
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
        layerData_IzzMR2 = 0.0;
        
        if (layerData.get() != 0) {
            CubicChebyshevMassDistribution::CoefficientType md_lD_coeff;
            CubicChebyshevMassDistribution::resize(md_lD_coeff,0);
            md_lD_coeff[0][0][0] = 0;
            osg::ref_ptr<CubicChebyshevMassDistribution> md_lD =
                new CubicChebyshevMassDistribution(md_lD_coeff,
                                                   R0_plate,
                                                   layerData.get());
            orsa::Cache<orsa::Vector> CM = sampled_CM;
            CM.lock();
            CCMD2SH(CM,
                    layerData_norm_C,
                    layerData_norm_S,
                    layerData_IzzMR2,
                    SH_degree, // gravityData->degree,
                    si.get(),
                    md_lD,
                    R0_plate,
                    R0_gravity);
            
            // scale for mass fraction later...
            
        }
    }
    
    double totalExcessMass=0.0;
    for (size_t k=0; k<ellipsoidLayerData.size(); ++k) {
        totalExcessMass += ellipsoidLayerData[k].excessMass;
    }
    const double layersTotalMassFraction = totalExcessMass / totalBodyMass;
    
    orsa::Cache<double> retVal;
    
    {
        const size_t i = row;
        
        orsa::Cache<double> layer_coeff;
        const QString ref_key = mod_gravityData_key(gravityData.get(),i);
        for (size_t l=0; l<=SH_degree; ++l) {
            for (size_t m=0; m<=l; ++m) {
                if (orsaPDS::RadioScienceGravityData::keyC(l,m) == ref_key) {
                    // ORSA_DEBUG("found: [%s] for l=%i, m=%i",ref_key.toStdString().c_str(),l,m);
                    layer_coeff        =     layersTotalMassFraction*layerData_norm_C[l][m].get_d();
                    // uniformShape_coeff = uniformShapeMassFraction*uniformShape_norm_C[l][m].get_d();
                } else if (orsaPDS::RadioScienceGravityData::keyS(l,m) == ref_key) {
                    // ORSA_DEBUG("found: [%s] for l=%i, m=%i",ref_key.toStdString().c_str(),l,m);
                    layer_coeff        =     layersTotalMassFraction*layerData_norm_S[l][m].get_d();
                    // uniformShape_coeff = uniformShapeMassFraction*uniformShape_norm_S[l][m].get_d();
                }                       
            }
        }
        if ("IzzMR2" == ref_key) {
            layer_coeff = layersTotalMassFraction*layerData_IzzMR2.get_d();
        }
        if (!layer_coeff.isSet()) {
            ORSA_DEBUG("problems...");
            exit(0);
        }
        
        // only for i == row
        retVal = layer_coeff;        
    }
    
    if (!retVal.isSet()) {
        ORSA_DEBUG("problems...");
        exit(0);
    }
    
    // ORSA_DEBUG("retVal: %g",(*retVal));
    
    return retVal;
}

