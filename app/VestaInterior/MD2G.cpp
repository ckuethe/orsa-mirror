#include <orsa/massDistribution.h>
#include <orsa/matrix.h>
#include <orsa/paulMoment.h>
#include <orsa/util.h>
#include <orsa/vector.h>

#include "vesta.h"

using namespace orsa;

int main() {
    
    osg::ref_ptr<VestaShape> shape = new VestaShape;
    if (!shape->read("vesta_thomas.dat")) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    osg::ref_ptr<orsa::MassDistribution> massDistribution = new orsa::UniformMassDistribution;

    orsa::Vector centerOfMass;
    orsa::Matrix shapeToLocal;
    orsa::Matrix localToShape;
    orsa::Matrix inertiaMatrix;
    osg::ref_ptr<orsa::PaulMoment> paulMoment;
          
    const unsigned int order = 4;
    const unsigned int N = 100000;

    const bool storeRandomVectors = true;
    
    osg::ref_ptr<RandomPointsInShape> randomPointsInShape =
        new RandomPointsInShape(shape,massDistribution,N,storeRandomVectors);
    
    const double volume = orsa::volume(randomPointsInShape);
    
    centerOfMass = orsa::centerOfMass(randomPointsInShape);
    // massDistribution);
    
    orsa::diagonalizedInertiaMatrix(shapeToLocal,
                                    localToShape,
                                    inertiaMatrix,
                                    centerOfMass,
                                    randomPointsInShape);
    // massDistribution);
    
    paulMoment = orsa::computePaulMoment(order,
                                         shapeToLocal,
                                         localToShape,
                                         centerOfMass,
                                         randomPointsInShape);
    // massDistribution);
    
	std::vector< std::vector<double> > C, S, norm_C, norm_S;
	std::vector<double> J;
	orsa::convert(C, S, norm_C, norm_S, J,
                  paulMoment, 
                  FromUnits(300,orsa::Unit::KM));
    
	ORSA_DEBUG("\\hline");
	ORSA_DEBUG("$x_{0}$    & $%+9.3f$ \\\\",orsa::FromUnits(centerOfMass.getX(),orsa::Unit::KM,-1));
	ORSA_DEBUG("$y_{0}$    & $%+9.3f$ \\\\",orsa::FromUnits(centerOfMass.getY(),orsa::Unit::KM,-1));
	ORSA_DEBUG("$z_{0}$    & $%+9.3f$ \\\\",orsa::FromUnits(centerOfMass.getZ(),orsa::Unit::KM,-1));
	ORSA_DEBUG("\%\\hline");
	for (unsigned int l=2; l<=order; ++l) {
        // J_l is minus C_l0, where C_l0 is not normalized
        ORSA_DEBUG("$J_{%i}$    & $%+9.6f$ \\\\",l,-C[l][0]);
	}
	ORSA_DEBUG("\%\\hline");
	for (unsigned int l=2; l<=order; ++l) {
        for (unsigned int m=0; m<=l; ++m) {
            // LaTeX Tabular style
            ORSA_DEBUG("$C_{%i%i}$   & $%+9.6f$ \\\\",l,m,norm_C[l][m]);
            if (m!=0) {
                ORSA_DEBUG("$S_{%i%i}$   & $%+9.6f$ \\\\",l,m,norm_S[l][m]);
            }
        }
	}
    
    
    return 0;
}
