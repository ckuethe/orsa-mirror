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

// DP = density profile

int main(int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if (argc != 6) {
        printf("Usage: %s <plate-model-file> <R0_km> <CCMDF-file> <num-sample-points> <equatorial_radius_km>\n",argv[0]);
        exit(0);
    }
    
    const double km   = orsa::FromUnits(1.0,orsa::Unit::KM);
    const double gcm3 = orsa::FromUnits(orsa::FromUnits(1.0,orsa::Unit::GRAM),orsa::Unit::CM,-3);
        
    const std::string plateModelFile = argv[1];
    const double plateModelR0 = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const std::string CCMDF_filename = argv[3];
    const int numSamplePoints = atoi(argv[4]);
    const double equatorialRadius = orsa::FromUnits(atof(argv[5]),orsa::Unit::KM);
    
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
    
    // rv = random vectors inside shape
    std::vector<orsa::Vector> rv;
    {
        const bool storeSamplePoints = false; // saving the points in rv
        osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape =
            new orsa::RandomPointsInShape(shapeModel,
        0,
        numSamplePoints,
        storeSamplePoints);
        orsa::Vector v;
        randomPointsInShape->reset();
        while (randomPointsInShape->get(v)) {
            rv.push_back(v);
        }
    }
    
    CubicChebyshevMassDistributionFile::DataContainer CCMDF;
    CubicChebyshevMassDistributionFile::read(CCMDF,CCMDF_filename);
    if (CCMDF.size() == 0) {
        ORSA_DEBUG("problem: empty CCMDF file");
        exit(0);
    }
    
    FILE * fp_DP_out = fopen("DP.out","w");
    FILE * fp_DP_surface_out = fopen("DP_surface.out","w");
    
    for (size_t ff=0; ff<CCMDF.size(); ++ff) {
        
        osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution = CCMD(CCMDF[ff]);
        
        std::vector< std::vector< std::vector<double> > > N;
        CCMD2ijk(N,
        2, // gravityData->degree,
        si.get(),
        massDistribution,
        plateModelR0);        
        
        // Izz / M R^2 using plate R
        const double inertiaMomentZZ_over_plateModelR0squared = N[2][0][0]/N[0][0][0] + N[0][2][0]/N[0][0][0] - N[1][0][0]/N[0][0][0] * N[1][0][0]/N[0][0][0] - N[0][1][0]/N[0][0][0] * N[0][1][0]/N[0][0][0];
        // alias
        const double CMa2 = inertiaMomentZZ_over_plateModelR0squared*orsa::square(plateModelR0)/orsa::square(equatorialRadius);
        
        {
            // density values at rv, to be sorted...
            std::list<double> ld;
            for (size_t k=0; k<rv.size(); ++k) {
                ld.push_back(massDistribution->density(rv[k]));
            }
                
            // sort
            ld.sort();
                
            double progress=1.0;
            orsa::Cache<double> old_progress, old_density;
            const double step=-1.0/ld.size();
            std::list<double>::const_iterator it_ld = ld.begin();
            while (it_ld != ld.end()) {
            
                const double density = (*it_ld);
            
                bool print=false;
                if (old_progress.isSet()) {
                    if ((cbrt(old_progress)-cbrt(progress))>0.001) {
                        print=true;
                    }
                } else {
                    print=true;
                }
                if (old_density.isSet()) {
                    if ((density-old_density)>0.01*gcm3) {
                        print=true;
                    }
                } else {
                    print=true;
                }
            
                if (print) {    
                    // volumetric version
                    // printf(fp_DP_out,"%.6f %.6f %.3f\n",CMa2,progress,density/gcm3);
                
                    // radial version
                    fprintf(fp_DP_out,"%.6f %.6f %.3f\n",CMa2,cbrt(progress),density/gcm3);
                
                    old_progress=progress;
                    old_density=density;
                }
            
                progress += step;
                ++it_ld;
                
            }
            fflush(fp_DP_out);
        }
    
        {
            osg::ref_ptr< orsa::Statistic<double> > stat_sv = new orsa::Statistic<double>;
            std::list<double> sd;
            for (size_t k=0; k<sv.size(); ++k) {
                const double density = massDistribution->density(sv[k]);
                stat_sv->insert(density);
                sd.push_back(density);
            }
            sd.sort();
            const double sd_size = sd.size();
            //
            double sd50_min, sd50_max;
            double sd90_min, sd90_max;
            double progress=0.0;
            double old_progress=progress;
            const double step=+1.0/sv.size();
            std::list<double>::const_iterator it_sd = sd.begin();
            while (it_sd != sd.end()) {
                const double density = (*it_sd);
                if (old_progress<=0.05 && progress>=0.05) {
                    sd90_min=density;
                }
                if (old_progress<=0.95 && progress>=0.95) {
                    sd90_max=density;
                }
                if (old_progress<=0.25 && progress>=0.25) {
                    sd50_min=density;
                }
                if (old_progress<=0.75 && progress>=0.75) {
                    sd50_max=density;
                }
                
                old_progress = progress;
                progress += step;
                ++it_sd;
            }
            //
            fprintf(fp_DP_surface_out,"%.6f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",
                CMa2,
                stat_sv->min()/gcm3,
                stat_sv->max()/gcm3,
                stat_sv->average()/gcm3,
                stat_sv->standardDeviation()/gcm3,
                sd50_min/gcm3,
                sd50_max/gcm3,
                sd90_min/gcm3,
                sd90_max/gcm3);
            fflush(fp_DP_surface_out);
        }
    }
    
    fclose(fp_DP_out);
    fclose(fp_DP_surface_out);
    
    return 0;
}
