#include "CubicChebyshevMassDistribution.h"
#include "simplex.h"
#include "gaskell.h"
#include <orsaPDS/RadioScienceGravity.h>

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
    
    CubicChebyshevMassDistributionFile::DataContainer CCMDF;
    CubicChebyshevMassDistributionFile::read(CCMDF,CCMDF_filename);
    
    osg::ref_ptr<orsaPDS::RadioScienceGravityData> gravityData = new orsaPDS::RadioScienceGravityData;
    orsaPDS::RadioScienceGravityFile::read(gravityData.get(),radioScienceGravityFile,512,1518);
    
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
            
#warning add output of inertia moments...
            
            gmp_fprintf(fp,
                        "%i %i %12.6f %12.6f %+12.6f %+12.6f %+12.6f %+12.6f %12.6f %8.6f %8.6f\n",
                        k,
                        elk,
                        flattening,
                        xy_flattening,
                        orsa::radToDeg()*lat_short_axis_pole,
                        orsa::radToDeg()*lon_short_axis_pole,
                        orsa::radToDeg()*lat_long_axis_pole,
                        orsa::radToDeg()*lon_long_axis_pole,
                        el->excessDensity,
                        el->volume()*el->excessDensity/(gravityData->GM/orsa::Unit::G()),
                        el->volume()/shapeModel->volume());
            fflush(fp);
            
            
        }
        
        
    }
    
    fclose(fp);
    
    return 0;
}
