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

int main(int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    /*
    if (argc != 7) {
        printf("Usage: %s <plate-model-file> <R0_km> <CCMDF-file> <rotation-period-hours> <num-sample-points> <epsAbs-km>\n",argv[0]);
        exit(0);
    }
    */
    //
    if (argc != 9) {
        printf("Usage: %s <plate-model-file> <R0_km> <CCMDF-file> <rotation-period-hours> <axis-tilt-deg> <axis-azimuth-deg> <num-sample-points> <epsAbs-km>\n",argv[0]);
        exit(0);
    }
    
    const double km   = orsa::FromUnits(1.0,orsa::Unit::KM);
    const double gcm3 = orsa::FromUnits(orsa::FromUnits(1.0,orsa::Unit::GRAM),orsa::Unit::CM,-3);
        
    const std::string plateModelFile = argv[1];
    const double plateModelR0 = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const std::string CCMDF_filename = argv[3];
    const double rotation_period = orsa::FromUnits(atof(argv[4]),orsa::Unit::HOUR);
    const double rotation_tilt   = orsa::degToRad()*atof(argv[5]);
    const double rotation_azim   = orsa::degToRad()*atof(argv[6]);
    const int numSamplePoints = atoi(argv[7]);
    const double epsAbs = orsa::FromUnits(atof(argv[8]),orsa::Unit::KM);
    
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
    
    CubicChebyshevMassDistributionFile::DataContainer CCMDF;
    CubicChebyshevMassDistributionFile::read(CCMDF,CCMDF_filename);
    if (CCMDF.size() == 0) {
        ORSA_DEBUG("problem: empty CCMDF file");
        exit(0);
    }
    //
    /*
    CubicChebyshevMassDistribution::CoefficientType coeff;
    CubicChebyshevMassDistribution::resize(coeff,x->T_degree); 
    //
    for (unsigned int i=0; i<=x->T_degree; ++i) {
        for (unsigned int j=0; j<=x->T_degree-i; ++j) {
            for (unsigned int k=0; k<=x->T_degree-i-j; ++k) {
                if (i+j+k<=x->T_degree) {
                    const size_t index = CubicChebyshevMassDistribution::index(i,j,k);
                    coeff[i][j][k] = gsl_vector_get(cT,index);
                }
            }            
        }
    }
    */
    //
    // CubicChebyshevMassDistributionFile::DataContainer::const_iterator it_CCMDF = CCMDF.end(); --it_CCMDF; // use the last one only...
    // osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution = CCMD(*it_CCMDF);
    osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution = CCMD(CCMDF[CCMDF.size()-1]);
     
    // sv = shape vectors = vectors on surface
    const std::vector<orsa::Vector> & sv = shapeModel->getVertexVector();
    
    // density value at sv
    std::vector<double> ds;
    ds.resize(sv.size());
    for (size_t k=0; k<sv.size(); ++k) {
        ds[k] = massDistribution->density(sv[k]);
    }
    
    {
        FILE * fp_ds = fopen("surface_density.xyz","w");
        ORSA_DEBUG("writing file [surface_density.xyz]...");
        for (size_t k=0; k<ds.size(); ++k) {
            const orsa::Vector & v = sv[k];
            const double lat = asin(v.getZ()/v.length());
            const double lon = fmod(orsa::twopi()+atan2(v.getY(),v.getX()),orsa::twopi());
            fprintf(fp_ds,"%g %g %g\n",
                    lon*orsa::radToDeg(),
                    lat*orsa::radToDeg(),
                    ds[k]/gcm3);
        }
        fclose(fp_ds);
    }
    
    if (0) {
        // density at random depth
        const size_t samples_per_vector = 1; // at least 1
        FILE * fp_dd = fopen("depth_density.dat","w");
        ORSA_DEBUG("writing file [depth_density.dat]...");
        for (size_t k=0; k<sv.size(); ++k) {
            const orsa::Vector p = sv[k]*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
            const double depth = (sv[k]-p).length();
            const double density = massDistribution->density(p);
            fprintf(fp_dd,"%g %g\n",depth/km,density/gcm3);
        }
        fclose(fp_dd);
    }
    
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
    
    // density value at rv
    std::vector<double> dv;
    dv.resize(rv.size());
    for (size_t k=0; k<rv.size(); ++k) {
        dv[k] = massDistribution->density(rv[k]);
    }
    {
        FILE * fp_dv = fopen("interior_density.dat","w");
        ORSA_DEBUG("writing file [interior_density.dat]...");
        for (size_t k=0; k<dv.size(); ++k) {
            fprintf(fp_dv,"%g\n",dv[k]/gcm3);
        }
        fclose(fp_dv);
    }
    
    // vR = virtualRadius
    const double vR  = cbrt(shapeModel->volume()/(rv.size()*(4.0*orsa::pi()/3.0)));
    const double vR2 = vR*vR;
    const double vR3 = vR*vR2;
    // const double OmegaSq = orsa::square(orsa::twopi()/rotation_period);
    // #warning need axis tilt info as input if nonzero...
    // const orsa::Vector omega = (orsa::twopi()/rotation_period)*orsa::Vector(0,0,1);
    const double tilt = rotation_tilt;
    const double azim = rotation_azim;
    const orsa::Vector omega = (orsa::twopi()/rotation_period)*orsa::Vector(sin(tilt)*cos(azim),sin(tilt)*sin(azim),cos(tilt));
    
    std::vector< std::vector< std::vector<double> > > N;
    CCMD2ijk(N,
             0, // degree
             si.get(),
             massDistribution,
             plateModelR0);
    const double totalMass = N[0][0][0]*orsa::cube(plateModelR0);
    ORSA_DEBUG("total mass: %g [kg]",orsa::FromUnits(totalMass,orsa::Unit::KG,-1));
    const double GM = orsa::Unit::G()*totalMass;
    
    // value of potential at sv; need this to compute average potential value to use as reference value for geoid
    std::vector<double> ps;
    ps.resize(sv.size());
    osg::ref_ptr< orsa::Statistic<double> > stat_ps = new orsa::Statistic<double>;
    for (size_t j=0; j<sv.size(); ++j) {
        ps[j] = massClusterPotential(sv[j],GM,omega,rv,dv,vR2,vR3);
        stat_ps->insert(ps[j]);
    }
    //
    const double pRef = stat_ps->average();
    ORSA_DEBUG("reference potential: %g [mks]",pRef);
    
    FILE * fp_geoid = fopen("geoid.xyz","w");
    ORSA_DEBUG("writing file [geoid.xyz]...");
    //
    FILE * fp_topo = fopen("topo_ref_geoid.xyz","w");
    ORSA_DEBUG("writing file [topo_ref_geoid.xyz]...");
    //
    // geoid vector
    std::vector<orsa::Vector> gv;
    gv.resize(sv.size());
    for (size_t j=0; j<sv.size(); ++j) {
        // search (bisection) for point with potential pRef along sv[j];
#warning could improve initial guess by using ps[j] above...
        const orsa::Vector u = sv[j].normalized();
        double lA = sv[j].length()*0.90;
        double pA = massClusterPotential(u*lA,GM,omega,rv,dv,vR2,vR3);
        double lB = sv[j].length()*1.10;
        double pB = massClusterPotential(u*lB,GM,omega,rv,dv,vR2,vR3);
        bool converged=false;
        if ((pRef>=std::min(pA,pB)) && (pRef<=std::max(pA,pB))) {
            while (1) {
                // ORSA_DEBUG("lA: %g pA: %g   lB: %g pB: %g",lA/km,pA,lB/km,pB);
                if (fabs(lA-lB)<epsAbs) { converged=true; break; }
                double lC = 0.5*(lA+lB);
                double pC = massClusterPotential(u*lC,GM,omega,rv,dv,vR2,vR3);
                if ((pRef>=std::min(pA,pC)) && (pRef<=std::max(pA,pC))) {
                    lB = lC;
                    pB = pC;
                } else {
                    lA = lC;
                    pA = pC;
                }
            }
        } else {
            ORSA_DEBUG("problems... probably you need to increase initial guessed range in code, or use more sampling points...");
            ORSA_DEBUG("lA: %g pA: %g   lB: %g pB: %g   pRef: %g",lA/km,pA,lB/km,pB,pRef);
            exit(0);
        }
        if (converged) {
            double lC = 0.5*(lA+lB);
            // double potential = massClusterPotential(u*lC,GM,OmegaSq,rv,dv,vR2,vR3);
            gv[j] = u*lC;
            
            // output
            const orsa::Vector & v = gv[j];
            const double lat = asin(v.getZ()/v.length());
            const double lon = fmod(orsa::twopi()+atan2(v.getY(),v.getX()),orsa::twopi());
            
            fprintf(fp_geoid,"%g %g %g\n",lon*orsa::radToDeg(),lat*orsa::radToDeg(),(gv[j].length())/km);
            fprintf(fp_topo, "%g %g %g\n",lon*orsa::radToDeg(),lat*orsa::radToDeg(),(sv[j].length()-gv[j].length())/km);
            fflush(fp_geoid);
            fflush(fp_topo);
        }
    }
    //
    fclose(fp_geoid);
    fclose(fp_topo);
       
    return 0;
}
