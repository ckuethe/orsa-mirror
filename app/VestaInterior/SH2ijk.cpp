#include "SH2ijk.h"

#include "vesta.h"
#include "gaskell.h"

#include "mpreal.h"

template <typename T> std::vector< std::vector< std::vector<size_t> > > SHIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SHIntegration<T>::index4Table;

int main(int argc, char **argv) {
    
    typedef mpfr::mpreal F;
    mpfr::mpreal::set_default_prec(256); // number of bits of precision
    
    /*
    if ( (argc != 7) &&
         (argc != 8) &&
         (argc != 9) ) {
        printf("Usage: %s <SH-model-file> <R0_km> <min-degree> <max-degree> <mod-N> <mod-i> [epsabs=0.0] [epsrel=0.0]\n",argv[0]);
        exit(0);
    }
    */
    if (argc != 7) {
        printf("Usage: %s <SH-model-file> <R0_km> <min-degree> <max-degree> <mod-N> <mod-i>\n",argv[0]);
        exit(0);
    }
    
    const std::string inputFile = argv[1];
    const double R0 = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const int minDegree = atoi(argv[3]);
    const int maxDegree = atoi(argv[4]);
    const int mod_N = atoi(argv[5]);
    const int mod_i = atoi(argv[6]);
    // const double epsabs = (argc >= 8) ? fabs(atof(argv[7])) : 0.0;
    // const double epsrel = (argc >= 9) ? fabs(atof(argv[8])) : 0.0;
    
    if ( (R0 <= 0.0) ||
         (minDegree < 0) ||
         (maxDegree < minDegree) ||
         (mod_N < 1) ||
         (mod_i < 0) ||
         (mod_i >= mod_N) ) {
        ORSA_DEBUG("invalid input...");
        exit(0);
    }

    // safer over NFS
    // sqlite3_vfs_register(sqlite3_vfs_find("unix-dotfile"), 1);
    
    const std::string SQLiteDBFileName = getSqliteDBFileName_SH(inputFile,R0);
    
    orsa::Debug::instance()->initTimer();
    
    // typedef std::vector< std::vector<double> > SHcoeff;
    SHcoeff norm_A;
    SHcoeff norm_B;
    readSH(norm_A,
           norm_B,
           inputFile,
           orsa::Unit::KM);
    
    // osg::ref_ptr<SHIntegration<F> > shi = new SHIntegration<F>(norm_A, norm_B, R0, epsabs, epsrel, SQLiteDBFileName);
    osg::ref_ptr<SHIntegration<F> > shi = new SHIntegration<F>(norm_A, norm_B, R0, SQLiteDBFileName);
    
    // shi->reserve(maxDegree);
    for (size_t degree=minDegree; degree<=(size_t)maxDegree; ++degree) {
        if ((degree%mod_N)==(size_t)mod_i) {
            for (size_t i=0; i<=degree; ++i) {
                for (size_t j=0; j<=degree; ++j) {
                    for (size_t k=0; k<=degree; ++k) {
                        if (i+j+k==degree) {
                            // const size_t index = SHIntegration<F>::getIndex(i,j,k);
                            const double integral_ijk = shi->getIntegral(i,j,k);
                            ORSA_DEBUG("integral [%02i,%02i,%02i]: %+20.12f",i,j,k,integral_ijk);
                        }
                    }
                }
            }
        }
    }
    
    return 0;
}
