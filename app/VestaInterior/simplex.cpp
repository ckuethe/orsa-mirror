#include "simplex.h"

#include "vesta.h"
#include "gaskell.h"

#include <qd/dd_real.h>
#include <qd/qd_real.h>

/*** CHOOSE ONE ***/
// typedef double T;
// typedef mpf_class T;
typedef dd_real T;
// typedef qd_real T;

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;

int main(int argc, char **argv) {
    
    if (argc != 7) {
        ORSA_DEBUG("Usage: %s <plate-model-file> <R0_km> <min-degree> <max-degree> <mod-N> <mod-i> ");
        exit(0);
    }
    
    const std::string inputFile = argv[1];
    const double R0 = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const int minDegree = atoi(argv[3]);
    const int maxDegree = atoi(argv[4]);
    const int mod_N = atoi(argv[5]);
    const int mod_i = atoi(argv[6]);
    
    if ( (R0 <= 0.0) ||
         (minDegree < 0) ||
         (maxDegree < minDegree) ||
         (mod_N < 1) ||
         (mod_i < 0) ||
         (mod_i >= mod_N) ) {
        ORSA_DEBUG("invalid input...");
        exit(0);
    }
    
    const std::string SQLiteDBFileName = getSqliteDBFileName(inputFile,R0);
    
    // QD
    unsigned int oldcw;
    fpu_fix_start(&oldcw);
    
    orsa::Debug::instance()->initTimer();
    
    //ORSA_DEBUG("current mpf precision: %i",mpf_get_default_prec());
    mpf_set_default_prec(128);
    // mpf_set_default_prec(256);
    // mpf_set_default_prec(512);
    // ORSA_DEBUG("updated mpf precision: %i",mpf_get_default_prec());
    
    /* osg::ref_ptr<VestaShape> vestaShape = new VestaShape;
       if (!vestaShape->read(inputFile)) {
       ORSA_ERROR("problems encountered while reading shape file...");
       exit(0);
       }
    */
    
    osg::ref_ptr<GaskellPlateModel> shapeModel = new GaskellPlateModel;
    if (!shapeModel->read(inputFile)) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    /* osg::ref_ptr<VestaShape> vestaShape = new VestaShape;
       if (!vestaShape->read("cube.dat")) {
       ORSA_ERROR("problems encountered while reading shape file...");
       exit(0);
       }
    */
    
    osg::ref_ptr<SimplexIntegration<T> > si = new SimplexIntegration<T>(shapeModel.get(), R0, SQLiteDBFileName);
    // osg::ref_ptr<SimplexIntegration> si_unit_R0 = new SimplexIntegration(vestaShape.get(),1.0);
    
    si->reserve(maxDegree);
    // si_unit_R0->reserve(maxDegree);
    for (size_t degree=minDegree; degree<=(size_t)maxDegree; ++degree) {
        for (size_t i=0; i<=degree; ++i) {
            for (size_t j=0; j<=degree; ++j) {
                for (size_t k=0; k<=degree; ++k) {
                    if (i+j+k==degree) {
                        const size_t index = SimplexIntegration<T>::getIndex(i,j,k);
                        if ((index%mod_N)==(size_t)mod_i) {

                            const double integral_ijk = si->getIntegral(i,j,k);
                            ORSA_DEBUG("integral [%02i,%02i,%02i]: %+20.12f",i,j,k,integral_ijk);
                            // ORSA_DEBUG("scaled: %+16.6e",integral_ijk*orsa::int_pow(R0,3+degree)); // 3=jacobian, degree=transformation
                        }
                    }
                }
            }
        }
    }
    
    // QD
    fpu_fix_end(&oldcw);
    
    return 0;
}
