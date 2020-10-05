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
#include "simplex.h"

#include "shape.h"

// #include "dislin.h"

using namespace std;
using namespace orsa;

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;

int main(int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if (argc != 4) {
        printf("Usage: %s <plate-model-file> <CCMDF-file> <step_km>\n",argv[0]);
        exit(0);
    }
    
    const std::string plateModelFile = argv[1];
    const std::string CCMDF_filename = argv[2];
    const double step_km = atof(argv[3]);
    
    // safer over NFS
    // sqlite3_vfs_register(sqlite3_vfs_find("unix-dotfile"), 1);
    
    osg::ref_ptr<InputShape> shapeModel = new InputShape;
       if (!shapeModel->read(plateModelFile)) {
       ORSA_ERROR("problems encountered while reading shape file...");
       exit(0);
    }
        
    CubicChebyshevMassDistributionFile::DataContainer CCMDF;
    CubicChebyshevMassDistributionFile::read(CCMDF,CCMDF_filename,1e12);
    if (CCMDF.size() == 0) exit(0);
    
    CubicChebyshevMassDistributionFile::DataContainer::const_iterator it_CCMDF = CCMDF.begin();
    
    // DUMP LAST ONE ONLY...
    it_CCMDF = CCMDF.end();
    --it_CCMDF;
    
    osg::ref_ptr<CubicChebyshevMassDistribution> md = CCMD(*it_CCMDF);
    
    char filename[1024];
    sprintf(filename,"dump_%gkm.3D",step_km);
    FILE * fp = fopen(filename,"w");
    fprintf(fp,"x y z density\n");

    /* 
       #warning USE THIS ONE!
       // USE THIS ONE!
       const double x_step = step_km;        
       const double x_min  = -350.00;
       const double x_max  =  350.00;
       //
       const double y_step = step_km;
       const double y_min  = -350.00;
       const double y_max  =  350.00;
       //
       const double z_step = step_km;
       const double z_min  = -350.00;
       const double z_max  =  350.00;
    */
    
#warning USE THIS ONE!
    // USE THIS ONE!
    const double x_step = step_km;        
    const double x_min  = -500.00;
    const double x_max  =  500.00;
    //
    const double y_step = step_km;
    const double y_min  = -500.00;
    const double y_max  =  500.00;
    //
    const double z_step = step_km;
    const double z_min  = -500.00;
    const double z_max  =  500.00;
    
    // 1D profile in X
    /* const double x_step = step_km;        
       const double x_min  = -80.00;
       const double x_max  = 100.00;
       //
       const double y_step = step_km;
       const double y_min  = -0.5*step_km;
       const double y_max  =  0.5*step_km;
       //
       const double z_step = step_km;
       const double z_min  = -0.5*step_km;
       const double z_max  =  0.5*step_km;
    */
    
    // 1D profile in Y
    /* const double x_step = step_km;
       const double x_min  = -0.5*step_km;
       const double x_max  =  0.5*step_km;
       //
       const double y_step = step_km;        
       const double y_min  = -80.00;
       const double y_max  = 100.00;
       //
       const double z_step = step_km;
       const double z_min  = -0.5*step_km;
       const double z_max  =  0.5*step_km;
    */
    
    // 1D profile in Z
    /* const double x_step = step_km;
       const double x_min  = -0.5*step_km;
       const double x_max  =  0.5*step_km;
       //
       const double y_step = step_km;
       const double y_min  = -0.5*step_km;
       const double y_max  =  0.5*step_km;
       //
       const double z_step = step_km;        
       const double z_min  = -80.00;
       const double z_max  = 100.00;
    */
    
    double x = x_min + 0.5*x_step;
    while (x<=x_max) {
        double y = y_min + 0.5*y_step;
        while (y<=y_max) {
            double z = z_min + 0.5*z_step;
            while (z<=z_max) {
                const orsa::Vector v(orsa::FromUnits(x,orsa::Unit::KM),
                                     orsa::FromUnits(y,orsa::Unit::KM),
                                     orsa::FromUnits(z,orsa::Unit::KM));
                
                if (shapeModel->isInside(v)) {
                    fprintf(fp,"%g %g %g %g\n",x,y,z,orsa::FromUnits(orsa::FromUnits(md->density(v),orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
                    fflush(fp);
                }
                
                z += z_step;
            }
            
            y += y_step;
        }
        
        x += x_step;
    }
    
    fclose(fp);
    
    ORSA_DEBUG("done.");
    
    return 0;
}
