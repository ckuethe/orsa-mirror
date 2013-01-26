#include "CubicChebyshevMassDistribution.h"
#include "simplex.h"
#include "gaskell.h"
#include <orsaPDS/RadioScienceGravity.h>
#include "CCMD2SH.h"
#include <orsa/statistic.h>

/*** CHOOSE ONE ***/
// typedef double simplex_T;
// typedef mpf_class simplex_T;
typedef dd_real simplex_T;
// typedef qd_real simplex_T;

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;

class index_value {
public:
    size_t index;
    double value;
public:
    static bool sort_by_value(const index_value & x, const index_value & y) {
        return x.value < y.value;
    }
};


// when using this index, IzzMR2 is used instead of one of the natural Clm or Slm in the gravityData file;
orsa::Cache<size_t> IzzMR2_index;
size_t get_IzzMR2_index() {
    if (IzzMR2_index.isSet()) {
        return IzzMR2_index;
    } else {
        ORSA_DEBUG("error: requesting unset value of IzzMR2_index");
        return -1;
    }
}
orsa::Cache<double> IzzMR2;
double get_IzzMR2() {
    if (IzzMR2.isSet()) {
        return IzzMR2;
    } else {
        ORSA_DEBUG("error: requesting unset value of IzzMR2");
        return 0.0;
    }
}
bool have_IzzMR2() {
    if (IzzMR2.isSet()) {
        return (IzzMR2 > 0.0);
    } else {
        return false;
    }
}

// modified versions of RadioScienceGravityData calls, to include C10,C11,S11 and IzzMR2
unsigned int mod_gravityData_index(const orsaPDS::RadioScienceGravityData * gravityData,
                                   const QString & key) {
    unsigned int index;
    if (key == orsaPDS::RadioScienceGravityData::keyC(0,0)) {
        index = 0;
    } else if (key == orsaPDS::RadioScienceGravityData::keyC(1,0)) {
        index = 1;
    } else if (key == orsaPDS::RadioScienceGravityData::keyC(1,1)) {
        index = 2;
    } else if (key == orsaPDS::RadioScienceGravityData::keyS(1,1)) {
        index = 3;
    } else if (key == "IzzMR2") {
        index = get_IzzMR2_index();
    } else {
        index = gravityData->index(key);
        if (index != 0) index += 3;
    }
    return index;
}

QString mod_gravityData_key(const orsaPDS::RadioScienceGravityData * gravityData,
                            const unsigned int & index) {
    // ORSA_DEBUG("index: %i   gravityData->numberOfCoefficients+3: %i",index,gravityData->numberOfCoefficients);
    if (index == 0) {
        return orsaPDS::RadioScienceGravityData::keyC(0,0);    
    } else if (index == 1) {
        return orsaPDS::RadioScienceGravityData::keyC(1,0);        
    } else if (index == 2) {
        return orsaPDS::RadioScienceGravityData::keyC(1,1);        
    } else if (index == 3) {
        return orsaPDS::RadioScienceGravityData::keyS(1,1);
    } else if (index == get_IzzMR2_index()) {
        return "IzzMR2";
    } else {
        return gravityData->key(index-3);
    }
}
//
double mod_gravityData_getCoeff(const orsaPDS::RadioScienceGravityData * gravityData,
                                const QString & key) {
    double coeff;
    if (key == orsaPDS::RadioScienceGravityData::keyC(0,0)) {
        coeff = 1.0;
    } else if ( (key == orsaPDS::RadioScienceGravityData::keyC(1,0)) ||  
                (key == orsaPDS::RadioScienceGravityData::keyC(1,1)) ||
                (key == orsaPDS::RadioScienceGravityData::keyS(1,1)) ) {
        coeff = 0.0;
    } else if (key == "IzzMR2") {
        coeff = get_IzzMR2();
    } else {
        coeff = gravityData->getCoeff(key);
    }
    return coeff;
}

double mod_gravityData_getCoeff(const orsaPDS::RadioScienceGravityData * gravityData,
                                const unsigned int & index) {
    return mod_gravityData_getCoeff(gravityData,mod_gravityData_key(gravityData,index));
}

unsigned int mod_gravityData_numberOfCoefficients(const orsaPDS::RadioScienceGravityData * gravityData) {
    // return (gravityData->numberOfCoefficients + 3 + (have_IzzMR2() ? 1 : 0));
    unsigned int retVal = gravityData->numberOfCoefficients + 3;
    if (get_IzzMR2_index()==retVal) ++retVal;
    return retVal;
}

gsl_vector * mod_gravityData_getCoefficientVector(const orsaPDS::RadioScienceGravityData * gravityData) {
    gsl_vector * mu = gravityData->getCoefficientVector();
    gsl_vector * mod_mu = gsl_vector_alloc(mod_gravityData_numberOfCoefficients(gravityData));
    for (unsigned int k=0; k<mod_gravityData_numberOfCoefficients(gravityData); ++k) {
        if (k==0) {
            // gsl_vector_set(mod_mu,k,gsl_vector_get(mu,k));
            gsl_vector_set(mod_mu,k,1.0);
        } else if ( (k==1) || (k==2) || (k==3) ) {
            gsl_vector_set(mod_mu,k,0.0);
        } else if (k == get_IzzMR2_index()) {
            gsl_vector_set(mod_mu,k,get_IzzMR2());
        } else {
            gsl_vector_set(mod_mu,k,gsl_vector_get(mu,k-3));
        }
    }
    gsl_vector_free(mu);
    return mod_mu;
}

gsl_matrix * mod_gravityData_getCovarianceMatrix(const orsaPDS::RadioScienceGravityData * gravityData) {
    gsl_matrix * covm = gravityData->getCovarianceMatrix();
    gsl_matrix * mod_covm = gsl_matrix_alloc(mod_gravityData_numberOfCoefficients(gravityData),
                                             mod_gravityData_numberOfCoefficients(gravityData));
    for (unsigned int l=0; l<mod_gravityData_numberOfCoefficients(gravityData); ++l) {
        for (unsigned int m=0; m<mod_gravityData_numberOfCoefficients(gravityData); ++m) { 
            if ((l==0) && (m==0)) {
                // gsl_matrix_set(mod_covm,l,m,gsl_matrix_get(covm,l,m));
                gsl_matrix_set(mod_covm,l,m,gsl_matrix_get(covm,l,m)/orsa::square(gravityData->GM));
#warning NEED TO SCALE THIS DOWN: from GM to 1.0 level... C_{00}... so check this!
            } else if ((l==1) || (l==2) || (l==3) || (l==get_IzzMR2_index()) || (m==1) || (m==2) || (m==3) || (m==get_IzzMR2_index())) {
                gsl_matrix_set(mod_covm,l,m,0.0);
            } else if ((l==0) && (m!=0)) {
                gsl_matrix_set(mod_covm,l,m,gsl_matrix_get(covm,l,m-3));
            } else if ((l!=0) && (m==0)) {
                gsl_matrix_set(mod_covm,l,m,gsl_matrix_get(covm,l-3,m));
            } else {
                gsl_matrix_set(mod_covm,l,m,gsl_matrix_get(covm,l-3,m-3));
            }
        }
    }
    gsl_matrix_free(covm);
    return mod_covm;   
}

gsl_matrix * mod_gravityData_getInverseCovarianceMatrix(const orsaPDS::RadioScienceGravityData * gravityData) {
    gsl_matrix * inv_covm = gravityData->getInverseCovarianceMatrix();
    gsl_matrix * mod_inv_covm = gsl_matrix_alloc(mod_gravityData_numberOfCoefficients(gravityData),
                                             mod_gravityData_numberOfCoefficients(gravityData));
    for (unsigned int l=0; l<mod_gravityData_numberOfCoefficients(gravityData); ++l) {
        for (unsigned int m=0; m<mod_gravityData_numberOfCoefficients(gravityData); ++m) { 
            if ((l==0) && (m==0)) {
                // gsl_matrix_set(mod_inv_covm,l,m,gsl_matrix_get(inv_covm,l,m));
                gsl_matrix_set(mod_inv_covm,l,m,gsl_matrix_get(inv_covm,l,m)/orsa::square(gravityData->GM));
#warning NEED TO SCALE THIS DOWN: from GM to 1.0 level... C_{00}... so check this!
            } else if ((l==1) || (l==2) || (l==3) || (l==get_IzzMR2_index()) || (m==1) || (m==2) || (m==3) || (m==get_IzzMR2_index())) {
                gsl_matrix_set(mod_inv_covm,l,m,0.0);
            } else if ((l==0) && (m!=0)) {
                gsl_matrix_set(mod_inv_covm,l,m,gsl_matrix_get(inv_covm,l,m-3));
            } else if ((l!=0) && (m==0)) {
                gsl_matrix_set(mod_inv_covm,l,m,gsl_matrix_get(inv_covm,l-3,m));
            } else {
                gsl_matrix_set(mod_inv_covm,l,m,gsl_matrix_get(inv_covm,l-3,m-3));
            }
        }
    }
    return mod_inv_covm;   
}

/**********/




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
    
    char line[4096];
    
    for (size_t k=0; k<CCMDF.size(); ++k) {
        
        // ORSA_DEBUG(" --- CCMDF %6i --- ",k);
        
        osg::ref_ptr<CubicChebyshevMassDistribution> massDistribution =
            new CubicChebyshevMassDistribution(CCMDF[k].coeff,
                                               plateModelR0,
                                               CCMDF[k].layerData);
        
        randomPointsInShape->updateMassDistribution(massDistribution.get());
        
        // some reference values; if any layer has different values, the whole solution is discarded
        orsa::Cache<double> ref_lat_short_axis_pole;
        
        bool discard=false;
        
        std::vector<index_value> ivv;
        for (size_t elk=0; elk<massDistribution->layerData->ellipsoidLayerVector.size(); ++elk) {
            index_value iv;
            iv.index = elk;
            iv.value = massDistribution->layerData->ellipsoidLayerVector[elk]->volume();
            ivv.push_back(iv);
        }
        std::sort(ivv.begin(),ivv.end(),index_value::sort_by_value);
        /* for (size_t elk=0; elk<massDistribution->layerData->ellipsoidLayerVector.size(); ++elk) {
           ORSA_DEBUG("ivv[%i] index: %i value: %f",elk,ivv[elk].index,ivv[elk].value);
           }
        */
        
        std::vector<double> total_density; // sum of all excess densities up to this layer
        std::vector<double> total_volume;  // layer volume minus the volume of the largest layer contained in this layer (it has a hole)
        std::vector<double> total_mass;    // mass of the layer, given by the product of the above total*density[j]*total_volume[j]
        double outer_density;
        double outer_mass;
        {
            outer_mass = gravityData->GM/orsa::Unit::G();
            for (size_t ivvj=0; ivvj<ivv.size(); ++ivvj) {
                outer_mass -=
                    massDistribution->layerData->ellipsoidLayerVector[ivv[ivvj].index]->excessDensity *
                    massDistribution->layerData->ellipsoidLayerVector[ivv[ivvj].index]->volume();
            }
            outer_density = outer_mass/shapeModel->volume();
            
#warning parameter here...
            if (outer_density < orsa::FromUnits(orsa::FromUnits(1.0,orsa::Unit::GRAM),orsa::Unit::CM,-3)) {
                ORSA_DEBUG("outer density too low, skipping...");
                discard=true;
                continue;
            }
            
            total_density.resize(ivv.size());
            for (size_t ivvj=0; ivvj<ivv.size(); ++ivvj) {
                total_density[ivv[ivvj].index] = outer_density;
                for (size_t ivvk=ivvj; ivvk<ivv.size(); ++ivvk) {
                    total_density[ivv[ivvj].index] += massDistribution->layerData->ellipsoidLayerVector[ivv[ivvk].index]->excessDensity;
                    // ORSA_DEBUG("ivvj: %i  ivvk: %i  total_density[%i]: %g",ivvj,ivvk,ivv[ivvj].index,total_density[ivv[ivvj].index]);
                }
            }
            total_volume.resize(ivv.size());
            for (size_t ivvj=0; ivvj<ivv.size(); ++ivvj) {
                total_volume[ivv[ivvj].index] = massDistribution->layerData->ellipsoidLayerVector[ivv[ivvj].index]->volume();
                if (ivvj>0) {
                    total_volume[ivv[ivvj].index] -= massDistribution->layerData->ellipsoidLayerVector[ivv[ivvj-1].index]->volume();
                }
                // ORSA_DEBUG("ivvj: %i  total_volume[%i]: %g    vol[%i] = %g",ivvj,ivv[ivvj].index,total_volume[ivv[ivvj].index],ivv[ivvj].index,massDistribution->layerData->ellipsoidLayerVector[ivv[ivvj].index]->volume());
            }
            total_mass.resize(ivv.size());
            for (size_t ivvj=0; ivvj<ivv.size(); ++ivvj) {
                total_mass[ivv[ivvj].index] = total_density[ivv[ivvj].index] * total_volume[ivv[ivvj].index];
            }
            /* for (size_t ivvj=0; ivvj<ivv.size(); ++ivvj) {
               ORSA_DEBUG("layer %i total_density: %g total_volume: %g total_mass: %g   vol[%i]: %g",
               ivv[ivvj].index,
               total_density[ivv[ivvj].index],
               total_volume[ivv[ivvj].index],
               total_mass[ivv[ivvj].index],
               ivv[ivvj].index,
               massDistribution->layerData->ellipsoidLayerVector[ivv[ivvj].index]->volume());
               }
            */
        }
        
        // init line
        
        gmp_sprintf(line,"%6i %10.3f ",
                    k,
                    orsa::FromUnits(orsa::FromUnits(outer_density,orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
        
        double old_a=0.0;
        double old_b=0.0;
        double old_c=0.0;
        for (size_t ivvk=0; ivvk<ivv.size(); ++ivvk) {
            const size_t elk = ivv[ivvk].index;
            // ORSA_DEBUG("elk: %i",elk);
            osg::ref_ptr<const LayerData::EllipsoidLayer> el = massDistribution->layerData->ellipsoidLayerVector[elk];
            const double a = el->a;
            const double b = el->b;
            const double c = el->c;
            
            if ( (a < old_a) ||
                 (b < old_b) ||
                 (c < old_c) ) {
                ORSA_DEBUG("intersecting layers, skipping...");
                discard=true;
                break;
            } else {
                old_a = a;
                old_b = b;
                old_c = c;
            }
            
            std::list<double> l;
            l.push_back(a);
            l.push_back(b);
            l.push_back(c);
            l.sort();
            std::list<double>::const_iterator it=l.begin();
            const double min_abc = (*it);
            ++it;
            const double mid_abc = (*it);
            ++it;
            const double max_abc = (*it);
            ++it;
            // ORSA_DEBUG("%g %g %g === %g %g %g",a,b,c,min_abc,mid_abc,max_abc);
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
            if (xy_flattening != 0.0) {
                ORSA_DEBUG("xy flattening is non-zero, skipping...");
                discard=true;
                break;
            }
            
            if (ref_lat_short_axis_pole.isSet()) {
                if (lat_short_axis_pole != ref_lat_short_axis_pole) {
                    ORSA_DEBUG("pole problems, skipping...");
                    discard=true;
                    break;
                }
            } else {
                ref_lat_short_axis_pole = lat_short_axis_pole;
            }
            
            char str[4096];
            gmp_sprintf(str,
                        "%2i %12.6f %12.6f %12.6f %+12.6f %+12.6f %+12.6f %+12.6f %+12.6f %+12.6f %.6f %.6f %.6f %12.6f %12.6f %+12.6f %+12.6f %+12.6f %+12.6f %10.3f %8.6f %8.6f %10.3f %8.6f %8.6f %12.6f ",
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
                        el->volume()/shapeModel->volume(),
                        orsa::FromUnits(orsa::FromUnits(total_density[elk],orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
                        total_mass[elk]/(gravityData->GM/orsa::Unit::G()),
                        total_volume[elk]/shapeModel->volume(),
                        orsa::FromUnits(cbrt(max_abc*mid_abc*min_abc),orsa::Unit::KM,-1));
            // sqrt(stat_D2->average()),
            // sqrt(stat_D4->average()),
            // sqrt(stat_D10->average()));
            strcat(line,str);
            
            /* ORSA_DEBUG("%g %g %g",
               el->volume(),
               total_volume[elk],
               shapeModel->volume());
            */
            
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
            // delta-gravity
            const size_t max_SH_degree = 4;
            std::vector< osg::ref_ptr< orsa::Statistic<double> > > stat;
            stat.resize(max_SH_degree+1);
            for (size_t shd=0; shd<stat.size(); ++shd) {
                stat[shd] = new orsa::Statistic<double>;
            }
            // osg::ref_ptr< orsa::Statistic<double> > stat_D2  = new orsa::Statistic<double>;
            // osg::ref_ptr< orsa::Statistic<double> > stat_D4  = new orsa::Statistic<double>;
            // osg::ref_ptr< orsa::Statistic<double> > stat_D10 = new orsa::Statistic<double>;
            {

                // need this just to compute CM
                std::vector< std::vector<mpf_class> > norm_C;
                std::vector< std::vector<mpf_class> > norm_S;
                orsa::Cache<orsa::Vector> CM;
                mpf_class IzzMR2;
                {
                    CM.reset();
                    CCMD2SH(CM,
                            norm_C,
                            norm_S,
                            IzzMR2,
                            2, // max_SH_degree, //gravityData->degree,
                            si.get(),
                            massDistribution.get(),
                            plateModelR0,
                            gravityData->R0);
                }
                CM.lock();
                
                // orsa::print(CM);
                
                std::vector< std::vector<mpf_class> > layerData_norm_C;
                std::vector< std::vector<mpf_class> > layerData_norm_S;
                mpf_class layerData_IzzMR2;
                {
                    CubicChebyshevMassDistribution::CoefficientType md_lD_coeff;
                    CubicChebyshevMassDistribution::resize(md_lD_coeff,0);
                    md_lD_coeff[0][0][0] = 0;
                    osg::ref_ptr<CubicChebyshevMassDistribution> md_lD =
                        new CubicChebyshevMassDistribution(md_lD_coeff,
                                                           plateModelR0,
                                                           massDistribution->layerData.get());
                    CM.lock();
                    CCMD2SH(CM,
                            layerData_norm_C,
                            layerData_norm_S,
                            layerData_IzzMR2,
                            max_SH_degree, // gravityData->degree,
                            si.get(),
                            md_lD,
                            plateModelR0,
                            gravityData->R0);
                    
                    const double layerMassFraction = massDistribution->layerData->totalExcessMass() / (gravityData->GM/orsa::Unit::G());
                    for (size_t l=0; l<=max_SH_degree; ++l) {
                        for (size_t m=0; m<=l; ++m) {
                            layerData_norm_C[l][m] *= layerMassFraction;
                            layerData_norm_S[l][m] *= layerMassFraction;
                        }
                    }
                    layerData_IzzMR2 *= layerMassFraction;
                }
                
                std::vector< std::vector<mpf_class> > NOLayerData_norm_C;
                std::vector< std::vector<mpf_class> > NOLayerData_norm_S;
                mpf_class NOLayerData_IzzMR2;
                {
#warning this is not just a no-layer coeff, it is a constant density coeff
                    CubicChebyshevMassDistribution::CoefficientType coeff;
                    CubicChebyshevMassDistribution::resize(coeff,0);
                    coeff[0][0][0] = 1.0;
                    osg::ref_ptr<CubicChebyshevMassDistribution> md_NOLD =
                        new CubicChebyshevMassDistribution(coeff,
                                                           plateModelR0,
                                                           0);
                    CM.lock();
                    CCMD2SH(CM,
                            NOLayerData_norm_C,
                            NOLayerData_norm_S,
                            NOLayerData_IzzMR2,
                            max_SH_degree, // SH_degree, // gravityData->degree,
                            si.get(),
                            md_NOLD.get(),
                        plateModelR0,
                            gravityData->R0);
                    const double NOLayerMassFraction = 1.0 - massDistribution->layerData->totalExcessMass() / (gravityData->GM/orsa::Unit::G());
                    for (size_t l=0; l<=max_SH_degree; ++l) {
                        for (size_t m=0; m<=l; ++m) {
                            NOLayerData_norm_C[l][m] *= NOLayerMassFraction;
                            NOLayerData_norm_S[l][m] *= NOLayerMassFraction;
                        }
                    }
                    NOLayerData_IzzMR2 *= NOLayerMassFraction;
                }
                
                for (size_t l=0; l<=max_SH_degree; ++l) {
                    for (size_t m=0; m<=l; ++m) {
                        const double gravity_data   = mod_gravityData_getCoeff(gravityData,orsaPDS::RadioScienceGravityData::keyC(l,m));
                        const double gravity_ccmd   = NOLayerData_norm_C[l][m].get_d();
                        const double gravity_layers = layerData_norm_C[l][m].get_d();
                        //
                        // const double delta = (gravity_data - gravity_ccmd - gravity_layers) / (gravity_data - gravity_ccmd);
                        const double delta = (gravity_data - gravity_ccmd - gravity_layers);
                        //
                        for (size_t shd=0; shd<stat.size(); ++shd) {
                            if (l<=shd) {
                                stat[shd]->insert(orsa::square(delta));
                            }
                        }
                        // if (l<=2)   stat_D2->insert(orsa::square(delta));
                        // if (l<=4)   stat_D4->insert(orsa::square(delta));
                        // if (l<=10) stat_D10->insert(orsa::square(delta));
                        // ORSA_DEBUG("C%i%i gravity_data: %g   gravity_ccmd: %g   gravity_layers: %g   delta: %g",l,m,gravity_data,gravity_ccmd,gravity_layers,delta);
                        if (m!=0) {
                            const double gravity_data   = mod_gravityData_getCoeff(gravityData,orsaPDS::RadioScienceGravityData::keyS(l,m));
                            const double gravity_ccmd   = NOLayerData_norm_S[l][m].get_d();
                            const double gravity_layers = layerData_norm_S[l][m].get_d();
                            //
                            // const double delta = (gravity_data - gravity_ccmd - gravity_layers) / (gravity_data - gravity_ccmd);
                            const double delta = (gravity_data - gravity_ccmd - gravity_layers);
                            //
                            for (size_t shd=0; shd<stat.size(); ++shd) {
                                if (l<=shd) {
                                    stat[shd]->insert(orsa::square(delta));
                                }
                            }
                            // if (l<=2)   stat_D2->insert(orsa::square(delta));
                            // if (l<=4)   stat_D4->insert(orsa::square(delta));
                            // if (l<=10) stat_D10->insert(orsa::square(delta));
                            // ORSA_DEBUG("S%i%i gravity_data: %g   gravity_ccmd: %g   gravity_layers: %g   delta: %g",l,m,gravity_data,gravity_ccmd,gravity_layers,delta);
                        }
                    }
                }                
            }
            
            char str[4096];
            for (size_t shd=0; shd<stat.size(); ++shd) {
                gmp_sprintf(str,
                            "%12.9f ",
                            sqrt(stat[shd]->sum()));
                strcat(line,str);
            }
        }        
        
        if (!discard) {
            gmp_fprintf(fp,"%s\n",line);
            fflush(fp);
            // ORSA_DEBUG("line: [%s]",line);
        }
        
    }
    
    fclose(fp);
    
    return 0;
}
