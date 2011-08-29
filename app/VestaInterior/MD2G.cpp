#include <orsa/massDistribution.h>
#include <orsa/matrix.h>
#include <orsa/paulMoment.h>
#include <orsa/util.h>
#include <orsa/vector.h>
#include <orsa/multifit.h>
#include <orsa/chebyshev.h> 

#include <orsaPDS/RadioScienceGravity.h>
#include "CubicChebyshevMassDistribution.h"
#include "simplex.h"

#include "vesta.h"
#include "gaskell.h"

#include <qd/dd_real.h>
#include <qd/qd_real.h>

using namespace orsa;

/*** CHOOSE ONE ***/
// typedef double simplex_T;
// typedef mpf_class simplex_T;
typedef dd_real simplex_T;
// typedef qd_real simplex_T;

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;


//// ChebyshevFit3D

const size_t ChebyshevFit3D_N = 4;

class ChebyshevFit3D : public orsa::Multifit {
public:
    ChebyshevFit3D() : 
        orsa::Multifit() { }
protected:
    void singleIterationDone(const orsa::MultifitParameters *) const {
        ORSA_DEBUG("--MARK--");
    }
protected:
    mutable std::vector< std::vector<double> >  Tx, Ty, Tz;
    void computeAllFunctionCalls(const orsa::MultifitParameters * /* par */, 
                                 const orsa::MultifitData       * data,
                                 const computeAllCallsMode        /* m */) const {
        // need to run this only once!
        static bool done=false;  
        if (!done) {
            const size_t N = ChebyshevFit3D_N; // copy from global var
            Tx.resize(data->size());
            Ty.resize(data->size());
            Tz.resize(data->size());
            for (size_t row=0; row<data->size(); ++row) {
                orsa::ChebyshevT(Tx[row],N,data->getD("x",row));
                orsa::ChebyshevT(Ty[row],N,data->getD("y",row));
                orsa::ChebyshevT(Tz[row],N,data->getD("z",row));
            }
            done=true;
        }
    }
public:
    double fun(const orsa::MultifitParameters * par, 
               const orsa::MultifitData       * /* data */,
               const unsigned int p, 
               const int          d,
               const unsigned int row) const {
        
        const size_t N = ChebyshevFit3D_N; // copy from global var
        
        double f = 0.0;
        
        // const double x = data->getD("x",row);

        // PRECOMPUTED!
        // std::vector<double> T;
        // orsa::ChebyshevT(T,N,x);
        
        char varName[1024];
        const size_t degree = N;
        for (size_t s=0; s<=degree; ++s) {
            for (size_t i=0; i<=degree; ++i) {
                for (size_t j=0; j<=degree-i; ++j) {
                    for (size_t k=0; k<=degree-i-j; ++k) {
                        if (i+j+k==s) {
                            sprintf(varName,"c%03i%03i%03i",i,j,k);
                            double c = par->get(varName);
                            if (p == par->index(varName)) {
                                c += d*par->getDelta(varName);
                            }
                            f += c*Tx[row][i]*Ty[row][j]*Tz[row][k];
                        }
                    }
                }
            }
        }
        
        return f;
    }
protected:
    int f_gsl(const gsl_vector * parameters, 
              void             * dataPoints, 
              gsl_vector       * f) {
    
        int retval =  orsa::Multifit::f_gsl(parameters, 
                                            dataPoints, 
                                            f);
        /* 
           const orsa::MultifitData * data = (orsa::MultifitData *) dataPoints;
           for (unsigned int j=0; j<data->size(); ++j) {
           ORSA_DEBUG("f[%02i] = %10.3f", j, gsl_vector_get(f,j));
           }
        */
        return retval;
    }
protected:
    int df_gsl (const gsl_vector * v, 
                void             * dataPoints, 
                gsl_matrix       * J) {
    
        int retval = Multifit::df_gsl(v, 
                                      dataPoints, 
                                      J);
    
        /* 
           const orsa::MultifitData * data = (orsa::MultifitData *) dataPoints;
           for (unsigned int k=0; k<_par->size(); ++k) {
           for (unsigned int j=0; j<data->size(); ++j) {
           ORSA_DEBUG("df[%02i][%02i] = %10.3f", j, k, gsl_matrix_get(J,j,k));
           }
           }
        */
        //
        return retval;
    }
  
};

void performChebyshevFit3D() {
    
    const size_t N = ChebyshevFit3D_N; // copy from global var
    
    osg::ref_ptr<orsa::MultifitParameters> par = 
        new orsa::MultifitParameters;
    char varName[1024];
    const size_t degree = N;
    for (size_t s=0; s<=degree; ++s) {
        for (size_t i=0; i<=degree; ++i) {
            for (size_t j=0; j<=degree-i; ++j) {
                for (size_t k=0; k<=degree-i-j; ++k) {
                    if (i+j+k==s) {
                        sprintf(varName,"c%03i%03i%03i",i,j,k);
                        par->insert(varName,0,1.0);
                    }
                }
            }
        }
    }
    
    // test ranges
    // par->setRange("a",40,42);
    // par->setRangeMin("a",40.22);
    // par->setRange("b",0,5);
    
    osg::ref_ptr<orsa::MultifitData> data = 
        new orsa::MultifitData;
    data->insertVariable("x");
    data->insertVariable("y");
    data->insertVariable("z");
    const double dx = 0.2; const size_t Nx = 2.0/dx;
    const double dy = 0.2; const size_t Ny = 2.0/dy;
    const double dz = 0.2; const size_t Nz = 2.0/dz;
    // const size_t maxRow = 2.0/dx;
    double f;
    // for (size_t row=0; row<=maxRow; ++row) {
    size_t row = 0, iter = 0;
    // while (1) {
    for (size_t kx=0; kx<=Nx; ++kx) {
        for (size_t ky=0; ky<=Ny; ++ky) {
            for (size_t kz=0; kz<=Nz; ++kz) {
                
                const double x = -1.0 + kx*dx;
                const double y = -1.0 + ky*dy;
                const double z = -1.0 + kz*dz;
                
                // ORSA_DEBUG("row: %i  iter: %i  x: %g y: %g z: %g",row,iter,x,y,z);
                
                // f = 1.0-0.3*x*x-4.0*x*y+z;
                f = 1.0+x+y+z+x*x+x*y+x*z+y*z+y*y+z*z;
                //
                data->insertD("x",row,x);
                data->insertD("y",row,y);
                data->insertD("z",row,z);
                data->insertF(row,f);
                data->insertSigma(row,1.0);
                ++row;
                ++iter;
            }
        }
    }

    const size_t maxRow = --row;
    char logFile[1024];
    snprintf(logFile,1024,"ChebyshevFit3D.log");
    
    osg::ref_ptr<ChebyshevFit3D> cf = new ChebyshevFit3D;
    //
    cf->setMultifitParameters(par.get());
    cf->setMultifitData(data.get());
    //
    cf->setLogFile(logFile);
    //
    cf->run();
    
    for (size_t row=0; row<=maxRow; ++row) {
        const double x = data->getD("x",row);
        const double y = data->getD("y",row);
        const double z = data->getD("z",row);
        const double f = data->getF(row);
        const double T = cf->fun(par.get(),
                                 data.get(),
                                 0,
                                 0,
                                 row);
        const double err = T-f;
        ORSA_DEBUG("FINAL: %g %g %g %g %g %g",
                   x,y,z,f,T,err);
    }
    
    for (size_t s=0; s<par->totalSize(); ++s) {
        ORSA_DEBUG("par[%03i] = [%s] = %+12.6f",s,par->name(s).c_str(),par->get(s));
    }
    
}

/************/


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
    
#warning add possibility to reduce the degree...
    
    if (argc != 5) {
        printf("Usage: %s <plate-model-file> <R0_km> <gravity-file-gravity-base-file> <output-gravity-file>\n",argv[0]);
        exit(0);
    }
    
    const std::string plateModelFile = argv[1];
    const double plateModelR0 = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const std::string radioScienceGravityBaseFile = argv[3];
    const std::string outputGravityFile = argv[4];
    
    if (plateModelR0 <= 0.0) {
        ORSA_DEBUG("invalid input...");
        exit(0);
    }
    
    // safer over NFS
    sqlite3_vfs_register(sqlite3_vfs_find("unix-dotfile"), 1);
    
    const std::string SQLiteDBFileName = getSqliteDBFileName(plateModelFile,plateModelR0);
    
    osg::ref_ptr<GaskellPlateModel> shapeModel = new GaskellPlateModel;
    if (!shapeModel->read(plateModelFile)) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    osg::ref_ptr<SimplexIntegration<simplex_T> > si =
        new SimplexIntegration<simplex_T>(shapeModel.get(), plateModelR0, SQLiteDBFileName);
    
    osg::ref_ptr<orsaPDS::RadioScienceGravityData> gravityData = new orsaPDS::RadioScienceGravityData;
    orsaPDS::RadioScienceGravityFile::read(gravityData.get(),radioScienceGravityBaseFile,512,1518);
    
    const double volume = si->getIntegral(0,0,0)*orsa::cube(plateModelR0);

    /* 
       const double CMx = si->getIntegral(1,0,0)*orsa::int_pow(plateModelR0,1);
       const double CMy = si->getIntegral(0,1,0)*orsa::int_pow(plateModelR0,1);
       const double CMz = si->getIntegral(0,0,1)*orsa::int_pow(plateModelR0,1);
       
       const double CMx_over_plateModelR0 = CMx / plateModelR0;
       const double CMy_over_plateModelR0 = CMy / plateModelR0;
       const double CMz_over_plateModelR0 = CMz / plateModelR0;
    */
    //
    const double alt_CMx_over_plateModelR0 = si->getIntegral(1,0,0) / si->getIntegral(0,0,0);
    const double alt_CMy_over_plateModelR0 = si->getIntegral(0,1,0) / si->getIntegral(0,0,0);
    const double alt_CMz_over_plateModelR0 = si->getIntegral(0,0,1) / si->getIntegral(0,0,0);
    
    // first determine the Chebyshev expansion of the mass distribution
    const size_t T_degree = 0;
    // const double gcm3 = orsa::FromUnits(orsa::FromUnits(1,orsa::Unit::GRAM),orsa::Unit::CM,-3);
    
#warning absolute or relative density? keep GM or change it?
    
    // (constant for now)
    CubicChebyshevMassDistribution::CoefficientType densityCCC; // CCC=CubicChebyshevCoefficient
    CubicChebyshevMassDistribution::resize(densityCCC,T_degree); 
    // densityCCC[0][0][0] = 3.4*gcm3;
    densityCCC[0][0][0] = 1.0;
    
    const double radiusCorrectionRatio = plateModelR0/gravityData->R0;
    
#warning track precision of operations
    
    for (size_t l=0; l<=gravityData->degree; ++l) {
        const double radiusCorrectionFactor = orsa::int_pow(radiusCorrectionRatio,l);
        for (size_t m=0; m<=l; ++m) {
            {
                const orsa::triIndex_mpq C_tri_integral = orsa::conversionCoefficients_C_integral(l,m);
                const orsa::triIndex_d   C_tri_norm     = orsa::conversionCoefficients_C_norm(l,m);
                double norm_C = 0.0;
                const orsa::triIndex_mpq S_tri_integral = orsa::conversionCoefficients_S_integral(l,m);
                const orsa::triIndex_d   S_tri_norm     = orsa::conversionCoefficients_S_norm(l,m);
                double norm_S = 0.0;
                // ni,nj,nk are the expansion of C_lm,S_lm in terms of N_ijk
                for (size_t ni=0; ni<=l; ++ni) {
                    for (size_t nj=0; nj<=l-ni; ++nj) {
                        for (size_t nk=0; nk<=l-ni-nj; ++nk) {
                            if ( (C_tri_integral[ni][nj][nk] == 0) && (S_tri_integral[ni][nj][nk] == 0) ) continue;
                            // ti,tj,tk are the expansion of the density in a cubic Chebyshev
                            for (size_t ti=0; ti<=T_degree; ++ti) {
                                for (size_t tj=0; tj<=T_degree-ti; ++tj) {
                                    for (size_t tk=0; tk<=T_degree-ti-tj; ++tk) {
                                        const std::vector<mpz_class> & cTi = orsa::ChebyshevTcoeff(ti);
                                        const std::vector<mpz_class> & cTj = orsa::ChebyshevTcoeff(tj);
                                        const std::vector<mpz_class> & cTk = orsa::ChebyshevTcoeff(tk);
                                        // ci,cj,ck are the expansion of each Chebyshev polinomial in terms of powers of x,y,z
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
                                                                    norm_C +=
                                                                        orsa::power_sign(bi+bj+bk) *
                                                                        densityCCC[ti][tj][tk] *
                                                                        C_tri_norm[ni][nj][nk] *
                                                                        mpz_class(orsa::binomial(ni,bi) *
                                                                                  orsa::binomial(nj,bj) *
                                                                                  orsa::binomial(nk,bk) *
                                                                                  cTi[ci] * cTj[cj] * cTk[ck]).get_d() *
                                                                        orsa::int_pow(alt_CMx_over_plateModelR0,bi) *
                                                                        orsa::int_pow(alt_CMy_over_plateModelR0,bj) *
                                                                        orsa::int_pow(alt_CMz_over_plateModelR0,bk) *
                                                                        si->getIntegral(ni-bi+ci,nj-bj+cj,nk-bk+ck);
                                                                }
                                                                
                                                                if (S_tri_integral[ni][nj][nk] != 0) {
                                                                    norm_S +=
                                                                        orsa::power_sign(bi+bj+bk) *
                                                                        densityCCC[ti][tj][tk] *
                                                                        S_tri_norm[ni][nj][nk] *
                                                                        mpz_class(orsa::binomial(ni,bi) *
                                                                                  orsa::binomial(nj,bj) *
                                                                                  orsa::binomial(nk,bk) *
                                                                                  cTi[ci] * cTj[cj] * cTk[ck]).get_d() *
                                                                        orsa::int_pow(alt_CMx_over_plateModelR0,bi) *
                                                                        orsa::int_pow(alt_CMy_over_plateModelR0,bj) *
                                                                        orsa::int_pow(alt_CMz_over_plateModelR0,bk) *
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
                        }
                    }
                }
                
                norm_C /= si->getIntegral(0,0,0);
                norm_C *= radiusCorrectionFactor;
                norm_S /= si->getIntegral(0,0,0);
                norm_S *= radiusCorrectionFactor;
                
                if (l>=2) {
                    gravityData->setCoeff(orsaPDS::RadioScienceGravityData::keyC(l,m),norm_C);
                    if (m!= 0) gravityData->setCoeff(orsaPDS::RadioScienceGravityData::keyS(l,m),norm_S);
                }
                
                ORSA_DEBUG("norm_C[%i][%i] = %g",l,m,norm_C);
                if (m != 0) ORSA_DEBUG("norm_S[%i][%i] = %g",l,m,norm_S);
            }
        }
    }
    
    ORSA_DEBUG("writing file [%s]",outputGravityFile.c_str());
    orsaPDS::RadioScienceGravityFile::write(gravityData.get(),outputGravityFile,512,1518);
    
    {
        // write barycenter file
        char filename[1024];
        sprintf(filename,"%s.barycenter_km.dat",outputGravityFile.c_str());
        FILE * fp = fopen(filename,"w");
        ORSA_DEBUG("writing file [%s]",filename);
        const double CMx = si->getIntegral(1,0,0)*plateModelR0;
        const double CMy = si->getIntegral(0,1,0)*plateModelR0;
        const double CMz = si->getIntegral(0,0,1)*plateModelR0;
        gmp_fprintf(fp,"%+9.3f %+9.3f %+9.3f\n",
                    orsa::FromUnits(CMx,orsa::Unit::KM,-1),
                    orsa::FromUnits(CMy,orsa::Unit::KM,-1),
                    orsa::FromUnits(CMz,orsa::Unit::KM,-1));
        fclose(fp);        
    }
    
    return 0;
}
