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
#include "penalty.h"
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

class ChebyshevFit3D : public orsa::Multifit {
public:
    ChebyshevFit3D(const size_t & T_degree) : 
        orsa::Multifit(),
        degree(T_degree) { }
protected:
    void singleIterationDone(const orsa::MultifitParameters *) const {
        ORSA_DEBUG("--MARK--");
    }
protected:
    const size_t degree;
    mutable std::vector< std::vector<double> >  Tx, Ty, Tz;
    void computeAllFunctionCalls(const orsa::MultifitParameters * /* par */, 
                                 const orsa::MultifitData       * data,
                                 const computeAllCallsMode        /* m */) const {
        // need to run this only once!
        static bool done=false;  
        if (!done) {
            Tx.resize(data->size());
            Ty.resize(data->size());
            Tz.resize(data->size());
            for (size_t row=0; row<data->size(); ++row) {
                orsa::ChebyshevT(Tx[row],degree,data->getD("x",row));
                orsa::ChebyshevT(Ty[row],degree,data->getD("y",row));
                orsa::ChebyshevT(Tz[row],degree,data->getD("z",row));
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
        
        double f = 0.0;
        
        // const double x = data->getD("x",row);
        
        // PRECOMPUTED!
        // std::vector<double> T;
        // orsa::ChebyshevT(T,N,x);
        
        char varName[1024];
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
    
#warning keep GM or change it?
    
    if ( (argc != 7) &&
         (argc != 8) ) {
        // passing CCMDF-input-file to use it as input mass distribution
        printf("Usage: %s <plate-model-file> <R0_km> <gravity-file-gravity-template-file> <output-gravity-file> <fitting-function-degree> <num-sample-points> [CCMDF-input-file]\n",argv[0]);
        exit(0);
    }
    
    const std::string plateModelFile = argv[1];
    const double plateModelR0 = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const std::string radioScienceGravityTemplateFile = argv[3];
    const std::string outputGravityFile = argv[4];
    const size_t T_degree_input = atoi(argv[5]);
    const size_t numSamplePoints = atoi(argv[6]);
    const bool have_CCMDF_file = (argc == 8);
    const std::string CCMDF_filename = (argc == 8) ? argv[7] : "";
    
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
    orsaPDS::RadioScienceGravityFile::read(gravityData.get(),radioScienceGravityTemplateFile,512,1518);
    
    const double GM = gravityData->GM; 
    const double volume = si->getIntegral(0,0,0)*orsa::cube(plateModelR0);
    const double bulkDensity = GM/orsa::Unit::G()/volume;
    
    CubicChebyshevMassDistribution::CoefficientType densityCCC; // CCC=CubicChebyshevCoefficient
    
    osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape;
    {
        const bool storeSamplePoints = true;
        randomPointsInShape =
            new orsa::RandomPointsInShape(shapeModel,
                                          0,
                                          numSamplePoints,
                                          storeSamplePoints);   
    }
    
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
        densityCCC = CCMDF[0].coeff;
        
    } else {
        
        // first determine the Chebyshev expansion of the mass distribution
        // #warning THIS MUST BE A PARAMETER
        const size_t T_degree = T_degree_input;
        
        // using relative density (coeff[0][0][0]=1 for constant density = bulk density)
        // CubicChebyshevMassDistribution::CoefficientType densityCCC; // CCC=CubicChebyshevCoefficient
        CubicChebyshevMassDistribution::resize(densityCCC,T_degree); 

        // choose mass distribution
        osg::ref_ptr<orsa::MassDistribution> massDistribution;
        //
        {
            const orsa::Vector coreCenter(orsa::FromUnits(0.0,orsa::Unit::KM),
                                          orsa::FromUnits(0.0,orsa::Unit::KM),
                                          orsa::FromUnits(0.0,orsa::Unit::KM));
            const double coreDensity = orsa::FromUnits(orsa::FromUnits(6.0,orsa::Unit::GRAM),orsa::Unit::CM,-3);
            // const double coreDensity = bulkDensity;
            const double mantleDensity = orsa::FromUnits(orsa::FromUnits(3.2,orsa::Unit::GRAM),orsa::Unit::CM,-3);
            // const double mantleDensity = bulkDensity;
            const double coreRadius = cbrt((3.0/(4.0*pi()))*volume*(bulkDensity-mantleDensity)/(coreDensity-mantleDensity));
            ORSA_DEBUG("coreRadius: %g [km]", orsa::FromUnits(coreRadius,orsa::Unit::KM,-1));
            massDistribution = new orsa::SphericalCorePlusMantleMassDistribution(coreCenter,
                                                                                 coreRadius,
                                                                                 coreDensity,
                                                                                 mantleDensity);
        }
        
        {
            // multifit of mass distribution
        
            osg::ref_ptr<orsa::MultifitParameters> par = 
                new orsa::MultifitParameters;
            char varName[1024];
            for (size_t s=0; s<=T_degree; ++s) {
                for (size_t i=0; i<=T_degree; ++i) {
                    for (size_t j=0; j<=T_degree-i; ++j) {
                        for (size_t k=0; k<=T_degree-i-j; ++k) {
                            if (i+j+k==s) {
                                sprintf(varName,"c%03i%03i%03i",i,j,k);
                                par->insert(varName,0,1.0);
                            }
                        }
                    }
                }
            }
            
            /* const bool storeSamplePoints = true;
               osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape =
               new orsa::RandomPointsInShape(shapeModel,
               massDistribution.get(),
               numSamplePoints,
               storeSamplePoints);       
            */
            //
            randomPointsInShape->updateMassDistribution(massDistribution.get());
            
            ORSA_DEBUG("MD penalty: %g",
                       MassDistributionPenalty(randomPointsInShape.get()));
            
            osg::ref_ptr<orsa::MultifitData> data = 
                new orsa::MultifitData;
            data->insertVariable("x");
            data->insertVariable("y");
            data->insertVariable("z");
            orsa::Vector v;
            double density;
            const double oneOverR0 = 1.0/plateModelR0;
            randomPointsInShape->reset();
            size_t row=0;
            size_t iter=0;
            while (randomPointsInShape->get(v,density)) { 
                data->insertD("x",row,v.getX()*oneOverR0);
                data->insertD("y",row,v.getY()*oneOverR0);
                data->insertD("z",row,v.getZ()*oneOverR0);
                data->insertF(row,density/bulkDensity);
                data->insertSigma(row,1.0);
                ++row;
                ++iter;
            }
            const size_t maxRow = --row;
        
            char logFile[1024];
            snprintf(logFile,1024,"MD2G_ChebyshevFit3D.log");
        
            osg::ref_ptr<ChebyshevFit3D> cf = new ChebyshevFit3D(T_degree);
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
                ORSA_DEBUG("FINAL: %+.3f %+.3f %+.3f %+.3f %+.3f %+.3f",
                           x,y,z,f,T,err);
            }
        
            for (size_t s=0; s<par->totalSize(); ++s) {
                ORSA_DEBUG("par[%03i] = [%s] = %+12.6f",s,par->name(s).c_str(),par->get(s));
            }
        
        
            for (size_t s=0; s<=T_degree; ++s) {
                for (size_t i=0; i<=T_degree; ++i) {
                    for (size_t j=0; j<=T_degree-i; ++j) {
                        for (size_t k=0; k<=T_degree-i-j; ++k) {
                            if (i+j+k==s) {
                                sprintf(varName,"c%03i%03i%03i",i,j,k);
                                densityCCC[i][j][k] = par->get(varName);
                            }
                        }
                    }
                }
            }
        }
    }
    
    const size_t T_degree = densityCCC.size()-1;
    
    orsa::Cache<double> penalty;
    {
        osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
            new CubicChebyshevMassDistribution(densityCCC,
                                               bulkDensity,     
                                               plateModelR0);
        randomPointsInShape->updateMassDistribution(massDistribution.get());
        penalty = MassDistributionPenalty(randomPointsInShape.get());
    }
    ORSA_DEBUG("CCC penalty: %g",(*penalty));
    
    {
        // another quick output...
#warning pass filename as parameter...
        CubicChebyshevMassDistributionFile::CCMDF_data data;
#warning review all these entries
        data.minDensity = 0.0;
        data.maxDensity = 0.0;
        data.deltaDensity = 0.0;
        data.penalty = 0.0;
        data.densityScale = bulkDensity;
        data.R0 = plateModelR0;
        data.SH_degree = gravityData->degree;
        data.coeff = densityCCC;
        CubicChebyshevMassDistributionFile::write(data,"MD2G.CCMDF.out");
    }
    
    const double radiusCorrectionRatio = plateModelR0/gravityData->R0;
    
    double i0d=0.0;
    double iXd=0.0;
    double iYd=0.0;
    double iZd=0.0;
    // ti,tj,tk are the expansion of the density in terms of the cubic Chebyshev 
    for (size_t ti=0; ti<=T_degree; ++ti) {
        for (size_t tj=0; tj<=T_degree-ti; ++tj) {
            for (size_t tk=0; tk<=T_degree-ti-tj; ++tk) {
                const std::vector<mpz_class> & cTi = orsa::ChebyshevTcoeff(ti);
                const std::vector<mpz_class> & cTj = orsa::ChebyshevTcoeff(tj);
                const std::vector<mpz_class> & cTk = orsa::ChebyshevTcoeff(tk);
                // ci,cj,ck are the expansion of each Chebyshev polynomial in terms of powers of x,y,z
                for (size_t ci=0; ci<=ti; ++ci) {
                    if (cTi[ci] == 0) continue;
                    for (size_t cj=0; cj<=tj; ++cj) {
                        if (cTj[cj] == 0) continue;
                        for (size_t ck=0; ck<=tk; ++ck) {
                            if (cTk[ck] == 0) continue;
                            const double baseFactor = densityCCC[ti][tj][tk] *
                                mpz_class(cTi[ci] * cTj[cj] * cTk[ck]).get_d();
                            i0d += baseFactor * si->getIntegral(ci,cj,ck);
                            iXd += baseFactor * si->getIntegral(ci+1,cj,ck);
                            iYd += baseFactor * si->getIntegral(ci,cj+1,ck);
                            iZd += baseFactor * si->getIntegral(ci,cj,ck+1);
                        }
                    }
                }
            }
        }
    }
    const double CMx_over_plateModelR0 = iXd / i0d;
    const double CMy_over_plateModelR0 = iYd / i0d;
    const double CMz_over_plateModelR0 = iZd / i0d;
    
    /* ORSA_DEBUG("si->getIntegral(0,0,0): %g",si->getIntegral(0,0,0));
       ORSA_DEBUG("i0d: %g",i0d);
       ORSA_DEBUG("iXd: %g",iXd);
       ORSA_DEBUG("iYd: %g",iYd);
       ORSA_DEBUG("iZd: %g",iZd);
    */
    
    {
        // write barycenter file
        char filename[1024];
        sprintf(filename,"%s.barycenter_km.dat",outputGravityFile.c_str());
        FILE * fp = fopen(filename,"w");
        ORSA_DEBUG("writing file [%s]",filename);
        const double CMx = CMx_over_plateModelR0*plateModelR0;
        const double CMy = CMy_over_plateModelR0*plateModelR0;
        const double CMz = CMz_over_plateModelR0*plateModelR0;
        gmp_fprintf(fp,"%+9.3f %+9.3f %+9.3f\n",
                    orsa::FromUnits(CMx,orsa::Unit::KM,-1),
                    orsa::FromUnits(CMy,orsa::Unit::KM,-1),
                    orsa::FromUnits(CMz,orsa::Unit::KM,-1));
        fclose(fp);        
    }
    
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
                            // ti,tj,tk are the expansion of the density in terms of the cubic Chebyshev
                            for (size_t ti=0; ti<=T_degree; ++ti) {
                                for (size_t tj=0; tj<=T_degree-ti; ++tj) {
                                    for (size_t tk=0; tk<=T_degree-ti-tj; ++tk) {
                                        const std::vector<mpz_class> & cTi = orsa::ChebyshevTcoeff(ti);
                                        const std::vector<mpz_class> & cTj = orsa::ChebyshevTcoeff(tj);
                                        const std::vector<mpz_class> & cTk = orsa::ChebyshevTcoeff(tk);
                                        // ci,cj,ck are the expansion of each Chebyshev polynomial in terms of powers of x,y,z
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
                                                                        orsa::int_pow(CMx_over_plateModelR0,bi) *
                                                                        orsa::int_pow(CMy_over_plateModelR0,bj) *
                                                                        orsa::int_pow(CMz_over_plateModelR0,bk) *
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
                        }
                    }
                }
                
                norm_C /= i0d;
                norm_C *= radiusCorrectionFactor;
                norm_S /= i0d;
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
    
    return 0;
}
