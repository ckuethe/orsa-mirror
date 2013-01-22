#include "CubicChebyshevMassDistribution.h"
#include "simplex.h"
#include "gaskell.h"
#include <orsaPDS/RadioScienceGravity.h>
#include "CCMD2SH.h"

/*** CHOOSE ONE ***/
// typedef double simplex_T;
// typedef mpf_class simplex_T;
typedef dd_real simplex_T;
// typedef qd_real simplex_T;

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;

int main(int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if (argc != 6) {
        printf("Usage: %s <RadioScienceGravityFile> <plate-model-file> <plate-model-R0_km> <num-sample-points> <CCMDF-input-file>\n",argv[0]);
        exit(0);
    }
    
    const std::string radioScienceGravityFile = argv[1];
    const std::string plateModelFile = argv[2];
    const double plateModelR0 = orsa::FromUnits(atof(argv[3]),orsa::Unit::KM);
    const int numSamplePoints = atoi(argv[4]);
    const std::string CCMDF_filename = argv[5];
    
    const std::string SQLiteDBFileName = getSqliteDBFileName_simplex(plateModelFile,plateModelR0);
    
    osg::ref_ptr<GaskellPlateModel> shapeModel = new GaskellPlateModel;
    if (!shapeModel->read(plateModelFile)) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    osg::ref_ptr<SimplexIntegration<simplex_T> > si = new SimplexIntegration<simplex_T>(shapeModel.get(), plateModelR0, SQLiteDBFileName);
    
    CubicChebyshevMassDistributionFile::DataContainer CCMDF;
    CubicChebyshevMassDistributionFile::read(CCMDF,CCMDF_filename);
    
    osg::ref_ptr<orsaPDS::RadioScienceGravityData> gravityData = new orsaPDS::RadioScienceGravityData;
    orsaPDS::RadioScienceGravityFile::read(gravityData.get(),radioScienceGravityFile,512,1518);
    
    if (0) {
        // test expected flatteing for a core at a given rotation period
        const double G = orsa::Unit::G();
        const double M = gravityData->GM/G;
        const double V = shapeModel->volume();
        const double T = orsa::FromUnits(5.342128,orsa::Unit::HOUR);
        const double Omega = orsa::twopi()/T;
        
        const double gcm3 = orsa::FromUnits(orsa::FromUnits(1.0,orsa::Unit::GRAM),orsa::Unit::CM,-3);
        
        double excess_mass_fraction = 0.0;
        while (excess_mass_fraction <= 1.0) {
            
            double excessDensity = 0.0;
            while (excessDensity <= 10.0*gcm3) {
                
                const double f =
                    15.0/(16.0*orsa::pi()) *
                    orsa::square(Omega)/G *
                    1.0/(excessDensity + (M/V)*(1.0-excess_mass_fraction));
                
                ORSA_DEBUG("EFEF %g %g %g %g",
                           T,
                           excess_mass_fraction,
                           excessDensity/gcm3,
                           f);
                
                excessDensity += 0.1*gcm3;
            }
            
            excess_mass_fraction += 0.01;
        }
    }
    
    osg::ref_ptr<orsa::RandomPointsInShape> randomPointsInShape;
    {
        const bool storeSamplePoints = true;
        randomPointsInShape =
            new orsa::RandomPointsInShape(shapeModel,
                                          0,
                                          numSamplePoints,
                                          storeSamplePoints);   
    }
    
    FILE * fp = fopen("scan.out","w");
    ORSA_DEBUG("writing file [scan.out]");
    
    for (size_t k=0; k<CCMDF.size(); ++k) {
        
        osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
            new CubicChebyshevMassDistribution(CCMDF[k].coeff,
                                               plateModelR0,
                                               CCMDF[k].layerData);
        
        randomPointsInShape->updateMassDistribution(massDistribution.get());
        
        char line[4096];
        
        gmp_sprintf(line,"%6i ",k);
        
        // some reference values; if any layer has different values, the whole solution is discarded
        orsa::Cache<double> ref_lat_short_axis_pole;
        
        bool discard=false;
        
        // flattening and short-axis pole orientation of ellipsoid layers
        for (size_t elk=0; elk<massDistribution->layerData->ellipsoidLayerVector.size(); ++elk) {
            osg::ref_ptr<const LayerData::EllipsoidLayer> el = massDistribution->layerData->ellipsoidLayerVector[elk];
            const double a = el->a;
            const double b = el->b;
            const double c = el->c;
            const double min_abc = std::min(a,std::min(b,c));
            const double max_abc = std::max(a,std::max(b,c));
            double tmp_mid=0.0;
            if ((a!=min_abc) && (a!=max_abc)) tmp_mid=a;
            if ((b!=min_abc) && (b!=max_abc)) tmp_mid=b;
            if ((c!=min_abc) && (c!=max_abc)) tmp_mid=c;
            if (tmp_mid==0.0) {
                ORSA_DEBUG("problems...");
                continue;
            }
            const double mid_abc = tmp_mid;
#warning review flattening definition: divide by max or by average?
            const double flattening = (max_abc-min_abc)/max_abc;
            const double xy_flattening = (max_abc-mid_abc)/max_abc;
            
            orsa::Vector u_tmp;
            
            if (a==min_abc) u_tmp = orsa::Vector(1,0,0);
            else if (b==min_abc) u_tmp = orsa::Vector(0,1,0);
            else if (c==min_abc) u_tmp = orsa::Vector(0,0,1);
            else {
                ORSA_DEBUG("problems...");
                continue;
            }
            //
            const orsa::Vector tmp_u_short_axis = el->rot * u_tmp;
            const orsa::Vector u_short_axis = (tmp_u_short_axis.getZ() > 0.0) ? tmp_u_short_axis : -tmp_u_short_axis;
            const double lat_short_axis_pole = orsa::halfpi()-acos(u_short_axis.getZ());
            const double lon_short_axis_pole = fmod(orsa::twopi()+atan2(u_short_axis.getY(),
                                                                        u_short_axis.getX()),orsa::twopi());
            
            if (a==max_abc) u_tmp = orsa::Vector(1,0,0);
            else if (b==max_abc) u_tmp = orsa::Vector(0,1,0);
            else if (c==max_abc) u_tmp = orsa::Vector(0,0,1);
            else {
                ORSA_DEBUG("problems...");
                continue;
            }
            const orsa::Vector tmp_u_long_axis = el->rot * u_tmp;
            const orsa::Vector u_long_axis = (tmp_u_long_axis.getZ() > 0.0) ? tmp_u_long_axis : -tmp_u_long_axis;
            const double lat_long_axis_pole = orsa::halfpi()-acos(u_long_axis.getZ());
            const double lon_long_axis_pole = fmod(orsa::twopi()+atan2(u_long_axis.getY(),
                                                                       u_long_axis.getX()),orsa::twopi());

            // inertia moments
            orsa::Cache<orsa::Vector> CM;
            mpf_class IxxMR2, IyyMR2, IzzMR2;
            inertia(CM,
                    IxxMR2,
                    IyyMR2,
                    IzzMR2,
                    // 2, // SH_degree,
                    si.get(),
                    massDistribution.get(),
                    plateModelR0);
            // orsa::print(CM);
            
            // filter?
            
            if (ref_lat_short_axis_pole.isSet()) {
                if (lat_short_axis_pole != ref_lat_short_axis_pole) {
                    discard=true;
                    break;
                }
            } else {
                ref_lat_short_axis_pole = lat_short_axis_pole;
            }
            
            /* gmp_fprintf(fp,
               "%6i %2i %12.6f %12.6f %12.6f %+12.6f %+12.6f %+12.6f %+12.6f %+12.6f %+12.6f %.6f %.6f %.6f %12.6f %12.6f %+12.6f %+12.6f %+12.6f %+12.6f %10.3f %8.6f %8.6f\n",
               k,
               elk,
               orsa::FromUnits(max_abc,orsa::Unit::KM,-1),
               orsa::FromUnits(mid_abc,orsa::Unit::KM,-1),
               orsa::FromUnits(min_abc,orsa::Unit::KM,-1),
               orsa::FromUnits(el->v0.getX(),orsa::Unit::KM,-1),
               orsa::FromUnits(el->v0.getY(),orsa::Unit::KM,-1),
               orsa::FromUnits(el->v0.getZ(),orsa::Unit::KM,-1),
               orsa::FromUnits((*CM).getX(),orsa::Unit::KM,-1),
               orsa::FromUnits((*CM).getY(),orsa::Unit::KM,-1),
               orsa::FromUnits((*CM).getZ(),orsa::Unit::KM,-1),
               IxxMR2.get_d(),
               IyyMR2.get_d(),
               IzzMR2.get_d(),
               flattening,
               xy_flattening,
               orsa::radToDeg()*lat_short_axis_pole,
               orsa::radToDeg()*lon_short_axis_pole,
               orsa::radToDeg()*lat_long_axis_pole,
               orsa::radToDeg()*lon_long_axis_pole,
               orsa::FromUnits(orsa::FromUnits(el->excessDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
               el->volume()*el->excessDensity/(gravityData->GM/orsa::Unit::G()),
               el->volume()/shapeModel->volume());
            */

            char str[4096];
            gmp_sprintf(str,
                        "%2i %12.6f %12.6f %12.6f %+12.6f %+12.6f %+12.6f %+12.6f %+12.6f %+12.6f %.6f %.6f %.6f %12.6f %12.6f %+12.6f %+12.6f %+12.6f %+12.6f %10.3f %8.6f %8.6f",
                        elk,
                        orsa::FromUnits(max_abc,orsa::Unit::KM,-1),
                        orsa::FromUnits(mid_abc,orsa::Unit::KM,-1),
                        orsa::FromUnits(min_abc,orsa::Unit::KM,-1),
                        orsa::FromUnits(el->v0.getX(),orsa::Unit::KM,-1),
                        orsa::FromUnits(el->v0.getY(),orsa::Unit::KM,-1),
                        orsa::FromUnits(el->v0.getZ(),orsa::Unit::KM,-1),
                        orsa::FromUnits((*CM).getX(),orsa::Unit::KM,-1),
                        orsa::FromUnits((*CM).getY(),orsa::Unit::KM,-1),
                        orsa::FromUnits((*CM).getZ(),orsa::Unit::KM,-1),
                        IxxMR2.get_d(),
                        IyyMR2.get_d(),
                        IzzMR2.get_d(),
                        flattening,
                        xy_flattening,
                        orsa::radToDeg()*lat_short_axis_pole,
                        orsa::radToDeg()*lon_short_axis_pole,
                        orsa::radToDeg()*lat_long_axis_pole,
                        orsa::radToDeg()*lon_long_axis_pole,
                        orsa::FromUnits(orsa::FromUnits(el->excessDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
                        el->volume()*el->excessDensity/(gravityData->GM/orsa::Unit::G()),
                        el->volume()/shapeModel->volume());
            strcat(line,str);
            
            if (0) {
                // test: output core-equatorial plane coordinates
                
                // core_excess_mass_fraction
                const double cf = el->volume()*el->excessDensity/(gravityData->GM/orsa::Unit::G());
                // selection conditions here...
                if (fabs(cf-0.25) < 0.001) {
                    
                    char filename[4096];
                    gmp_sprintf(filename,"core_equator_%8.6f.latlon",cf);
                    FILE * fp = fopen(filename,"w");
                    ORSA_DEBUG("writing file [%s]",filename);
                    const orsa::Vector u_X = u_long_axis;
                    const orsa::Vector u_Z = u_short_axis;
                    const orsa::Vector u_Y = orsa::externalProduct(u_Z,u_X);
                    double alpha = 0.0;
                    const double alpha_step = orsa::pi()/36.0;
                    double s,c;
                    orsa::Vector intersectionPoint;
                    orsa::Vector normal;
                    while (alpha <= orsa::twopi()+0.1*alpha_step) {
                        ::sincos(alpha,&s,&c);
                        orsa::Vector u = c*u_X + s*u_Y;
                        const bool good = shapeModel->rayIntersection(intersectionPoint,
                                                                      normal,
                                                                      el->v0,
                                                                      u,
                                                                      false);
                        /* if (!good) {
                           ORSA_DEBUG("problems...");
                           }
                        */
                        const orsa::Vector u_iP = intersectionPoint.normalized();
                        const double lat_iP = orsa::halfpi()-acos(u_iP.getZ());
                        const double lon_iP = fmod(orsa::twopi()+atan2(u_iP.getY(),
                                                                       u_iP.getX()),orsa::twopi());
                        
                        fprintf(fp,"%+12.6f %+12.6f\n",
                                orsa::radToDeg()*lat_iP,
                                orsa::radToDeg()*lon_iP);                        
                        alpha += alpha_step;
                    }
                    fclose(fp);
                    exit(0);
                }
                
            }

        }

        if (!discard) {
            gmp_fprintf(fp,"%s\n",line);
            fflush(fp);
        }
        
    }
    
    fclose(fp);
    
    return 0;
}
