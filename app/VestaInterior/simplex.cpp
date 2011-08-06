#include "simplex.h"

#include "vesta.h"

SimplexIntegration::IndexTableType  SimplexIntegration::indexTable;
SimplexIntegration::Index4TableType SimplexIntegration::index4Table;


int main() {
    
    orsa::Debug::instance()->initTimer();
    
    /* osg::ref_ptr<VestaShape> vestaShape = new VestaShape;
       if (!vestaShape->read("vesta_thomas.dat")) {
       ORSA_ERROR("problems encountered while reading shape file...");
       exit(0);
       }
    */
    
    osg::ref_ptr<VestaShape> vestaShape = new VestaShape;
    if (!vestaShape->read("cube.dat")) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    osg::ref_ptr<SimplexIntegration> si = new SimplexIntegration(vestaShape.get());


    const size_t maxDegree = 20;
    si->reserve(maxDegree);
    for (size_t degree=0; degree<=maxDegree; ++degree) {
        for (size_t i=0; i<=degree; ++i) {
            for (size_t j=0; j<=degree; ++j) {
                for (size_t k=0; k<=degree; ++k) {
                    if (i+j+k==degree) {
                        ORSA_DEBUG("integral [%02i,%02i,%02i]: %+16.6e",i,j,k,si->getIntegral(i,j,k));
                    }
                }
            }
        }
    }
    
    return 0;
}
