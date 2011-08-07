#include "simplex.h"

#include "vesta.h"

SimplexIntegration::IndexTableType  SimplexIntegration::indexTable;
SimplexIntegration::Index4TableType SimplexIntegration::index4Table;


int main() {
    
    orsa::Debug::instance()->initTimer();
    
    ORSA_DEBUG("current mpf precision: %i",mpf_get_default_prec());
    mpf_set_default_prec(128);
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
    
    osg::ref_ptr<SimplexIntegration> si = new SimplexIntegration(vestaShape.get(), R0);
    // osg::ref_ptr<SimplexIntegration> si_unit_R0 = new SimplexIntegration(vestaShape.get(),1.0);
    
    const size_t maxDegree = 200;
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
    
    return 0;
}
