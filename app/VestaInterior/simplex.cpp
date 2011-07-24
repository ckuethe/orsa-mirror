#include "simplex.h"

#include "vesta.h"

SimplexIntegration::IndexTableType SimplexIntegration::indexTable;

int main() {
    
    orsa::Debug::instance()->initTimer();
    
    osg::ref_ptr<VestaShape> vestaShape = new VestaShape;
    if (!vestaShape->read("vesta_thomas.dat")) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    osg::ref_ptr<SimplexIntegration> si = new SimplexIntegration(vestaShape.get());

    ORSA_DEBUG("integral [0,0,0]: %g",si->getIntegral(0,0,0));
    
    return 0;
}
