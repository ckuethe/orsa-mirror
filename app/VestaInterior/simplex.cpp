#include "simplex.h"

#include "vesta.h"

#include <qd/dd_real.h>
#include <qd/qd_real.h>

// Choose one
// typedef double T;
typedef mpf_class T;
// typedef dd_real T;
// typedef qd_real T;

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;

int main() {

    // QD
    unsigned int oldcw;
    fpu_fix_start(&oldcw);
    
    orsa::Debug::instance()->initTimer();
    
    ORSA_DEBUG("current mpf precision: %i",mpf_get_default_prec());
    // mpf_set_default_prec(128);
    // mpf_set_default_prec(256);
    // mpf_set_default_prec(512);
    ORSA_DEBUG("updated mpf precision: %i",mpf_get_default_prec());
    
    /* osg::ref_ptr<VestaShape> vestaShape = new VestaShape;
       if (!vestaShape->read("vesta_thomas.dat")) {
       ORSA_ERROR("problems encountered while reading shape file...");
       exit(0);
       }
       const double R0 = orsa::FromUnits(300.0,orsa::Unit::KM);
    */
    
    osg::ref_ptr<VestaShape> vestaShape = new VestaShape;
    if (!vestaShape->read("cube.dat")) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    const double R0 = orsa::FromUnits(1.0,orsa::Unit::KM);
    
    osg::ref_ptr<SimplexIntegration<T> > si = new SimplexIntegration<T>(vestaShape.get(), R0);
    // osg::ref_ptr<SimplexIntegration> si_unit_R0 = new SimplexIntegration(vestaShape.get(),1.0);
    
    const size_t maxDegree = 10;
    si->reserve(maxDegree);
    // si_unit_R0->reserve(maxDegree);
    for (size_t degree=0; degree<=maxDegree; ++degree) {
        for (size_t i=0; i<=degree; ++i) {
            for (size_t j=0; j<=degree; ++j) {
                for (size_t k=0; k<=degree; ++k) {
                    if (i+j+k==degree) {
                        // ORSA_DEBUG("integral [%02i,%02i,%02i]: %+16.6e",i,j,k,si->getIntegral(i,j,k));
                        ORSA_DEBUG("integral [%02i,%02i,%02i]: %+16.9g",i,j,k,si->getIntegral(i,j,k));
                        // ORSA_DEBUG("integral [%02i,%02i,%02i]: %+16.9e   [R0=1.0]",i,j,k,si_unit_R0->getIntegral(i,j,k));
                        // ORSA_DEBUG("scaled: %+16.6e",si->getIntegral(i,j,k)*orsa::int_pow(R0,3+degree)); // 3=jacobian, degree=transformation
                    }
                }
            }
        }
    }

    // QD
    fpu_fix_end(&oldcw);
    
    return 0;
}
