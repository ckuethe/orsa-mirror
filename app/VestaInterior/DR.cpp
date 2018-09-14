#include <orsa/massDistribution.h>
#include <orsa/matrix.h>
#include <orsa/paulMoment.h>
#include <orsa/util.h>
#include <orsa/vector.h>
#include <orsa/multifit.h>
#include <orsa/chebyshev.h> 
#include <orsa/debug.h>
#include <orsa/print.h>
#include <orsa/statistic.h>

#include <orsaPDS/RadioScienceGravity.h>
#include "CubicChebyshevMassDistribution.h"
#include "VestaInteriorAnalytic.h"
#include "simplex.h"
#include "CCMD2ijk.h"

#include "vesta.h"
#include "gaskell.h"
#include "eros_shape.h"

#include <libgen.h>

using namespace std;
using namespace orsa;

typedef mpfr::mpreal F;

// DR = depth range = range of depth with density within a given range

int main(int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if (argc != 6) {
        printf("Usage: %s <plate-model-file> <R0_km> <CCMDF-file> <rho_min_gcm3> <rho_max_gcm3>\n",argv[0]);
        exit(0);
    }
    
    const double km   = orsa::FromUnits(1.0,orsa::Unit::KM);
    const double gcm3 = orsa::FromUnits(orsa::FromUnits(1.0,orsa::Unit::GRAM),orsa::Unit::CM,-3);
        
    const std::string plateModelFile = argv[1];
    const double plateModelR0 = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const std::string CCMDF_filename = argv[3];
    const double rho_min = std::min(atof(argv[4]),atof(argv[5]))*gcm3;
    const double rho_max = std::max(atof(argv[4]),atof(argv[5]))*gcm3;
    
    if (plateModelR0 <= 0.0){
        ORSA_DEBUG("invalid input...");
        exit(0);
    }
    
    /* osg::ref_ptr<VestaShape> shapeModel = new VestaShape;
    if (!shapeModel->read(plateModelFile)) {
    ORSA_ERROR("problems encountered while reading shape file...");
    exit(0);
    }
    */
    
    /* 
    osg::ref_ptr<ErosShape> shapeModel = new ErosShape;
    if (!shapeModel->read(plateModelFile)) {
    ORSA_ERROR("problems encountered while reading shape file...");
    exit(0);
    }
    */
    
    osg::ref_ptr<GaskellPlateModel> shapeModel = new GaskellPlateModel;
    if (!shapeModel->read(plateModelFile)) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    const std::string SQLiteDBFileName = getSqliteDBFileName_simplex(plateModelFile,plateModelR0);
    osg::ref_ptr<SimplexIntegration<F> > si = new SimplexIntegration<F>(shapeModel.get(), plateModelR0, SQLiteDBFileName);    
    
    // sv = shape vectors = vectors on surface
    const std::vector<orsa::Vector> & sv = shapeModel->getVertexVector();
    
    CubicChebyshevMassDistributionFile::DataContainer CCMDF;
    CubicChebyshevMassDistributionFile::read(CCMDF,CCMDF_filename);
    if (CCMDF.size() == 0) {
        ORSA_DEBUG("problem: empty CCMDF file");
        exit(0);
    }
    if (CCMDF.size()!=1) {
        ORSA_DEBUG("only one entry per CCMDF...");
        exit(0);
    }
        
    osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution = CCMD(CCMDF[0]);
    
    // radial stepping size...
    const double ds = 0.1*km;
    
    for (size_t k=0; k<sv.size(); ++k) {
        const orsa::Vector r = sv[k];
        const double       l = r.length();
        const orsa::Vector u = r.normalized();
        osg::ref_ptr< orsa::Statistic<double> > stat_depth = new orsa::Statistic<double>;
        size_t j=0;
        while (j*ds<l) {
            const double depth = j*ds;
            const orsa::Vector p = r - depth*u;
            const double density = massDistribution->density(p);
            if ((density>=rho_min) && (density<=rho_max)) {
                stat_depth->insert(depth);
            }
            ++j;
        }
        // print
        const double lat = asin(u.getZ());
        const double lon = fmod(orsa::twopi()+atan2(u.getY(),u.getX()),orsa::twopi());
        printf("%7.3f %+7.3f %g %g %g\n",
            lon*orsa::radToDeg(),
            lat*orsa::radToDeg(),
            stat_depth->min()/km,
            stat_depth->average()/km,
            (stat_depth->max()-stat_depth->min())/km);
    }
    
    return 0;
}
