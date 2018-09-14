#include <orsa/debug.h>
#include <orsa/multimin.h>
#include <orsa/print.h>
#include <orsa/shape.h>
#include <orsa/statistic.h>
#include <orsa/util.h>
#include <orsa/vector.h>

class vertex {
public:
    orsa::Cache<orsa::Vector> v;
    orsa::Cache<double> lon, lat;
    orsa::Cache<bool> use;
public:
    void print() const {
        orsa::print(v);
        orsa::print(lon);
        orsa::print(lat);
        orsa::print(use);
    }
};

// #warning NOTE: should use the plumb line instead of the closestVertex, or maybe simply radial since it is in lat/lon
#warning REMEMBER TO USE the appropriate scaling when lat/lon points are not uniformly distributed

#warning can determine the uncertainty on fitting parameters?

class GeoidFitMultimin : public orsa::Multimin {
public:
    typedef std::vector<vertex> VertexVector;
public:
    GeoidFitMultimin(const VertexVector & vv) : orsa::Multimin(), vertexVector(vv) { }
public:
    const VertexVector vertexVector;
public:
    // both are in the same ref. sys. as the original vertexVector
    static void getEllipsoidVV(VertexVector & ellipsoidVertexVector,
                               const VertexVector & vertexVector,
                               const orsa::MultiminParameters * par,
                               const bool useAll=false) {
        osg::ref_ptr<orsa::EllipsoidShape> shape = new orsa::EllipsoidShape(par->get("a"),par->get("b"),par->get("c"));
        /*
        const orsa::Vector v2e_dv = orsa::Vector(par->get("v2e_dx"),par->get("v2e_dy"),par->get("v2e_dz"));
        orsa::Matrix m;
        if (!orsa::eulerAnglesToMatrix(m,
                                       par->get("v2e_psi"),
                                       par->get("v2e_theta"),
                                       par->get("v2e_phi"))) {
            ORSA_DEBUG("problems...");
            exit(0);
        }
        const orsa::Matrix v2e = m;
        if (!orsa::Matrix::invert(v2e,m)) {
            ORSA_DEBUG("problems...");
            exit(0);
        }
        const orsa::Matrix e2v = m;
        */
        ellipsoidVertexVector.resize(vertexVector.size());
        orsa::Vector intersectionPoint;
        orsa::Vector normal;
        for (size_t k=0; k<vertexVector.size(); ++k) {
            ellipsoidVertexVector[k].use = vertexVector[k].use;
            if (useAll || vertexVector[k].use) {
                
#warning CHOOSE one method here
                
                // first method: closestVertex
                // ellipsoidVertexVector[k].v = v2e_dv + e2v*(shape->closestVertex(v2e*(vertexVector[k].v-v2e_dv)));
                
                // second method: rayIntersection
                if (!shape->rayIntersection(intersectionPoint,
                                            normal,
                                            orsa::Vector(0,0,0),
                                            (*vertexVector[k].v).normalized(), // (v2e*(vertexVector[k].v-v2e_dv)).normalized(),
                                            false)) {
                    ORSA_DEBUG("problems...");
                    exit(0);
                }
                // ellipsoidVertexVector[k].v = v2e_dv + e2v*intersectionPoint;
                ellipsoidVertexVector[k].v = intersectionPoint;
                
                // derived quantities (don't comment out)
                ellipsoidVertexVector[k].lon = fmod(orsa::twopi()+atan2((*ellipsoidVertexVector[k].v).getY(),(*ellipsoidVertexVector[k].v).getX()),orsa::twopi());
                ellipsoidVertexVector[k].lat = asin((*ellipsoidVertexVector[k].v).normalized().getZ());
                
                // debug
                // vertexVector[k].print();
                // ellipsoidVertexVector[k].print();
                // ORSA_DEBUG("sp: %g",(*vertexVector[k].v).normalized()*(*ellipsoidVertexVector[k].v).normalized());
            }
        }
    }
public:
    double fun(const orsa::MultiminParameters * par) const {
        VertexVector ellipsoidVertexVector;
        getEllipsoidVV(ellipsoidVertexVector,
                       vertexVector,
                       par);
        double retVal = 0.0;
        size_t num = 0;
        for (size_t k=0; k<vertexVector.size(); ++k) {
            if (ellipsoidVertexVector[k].use) {
                retVal += (ellipsoidVertexVector[k].v-vertexVector[k].v).lengthSquared();
                ++num;
            }
        }
        retVal /= num;
        retVal = sqrt(retVal);

        /*** tmp: just to compute offset Z axis... ***/
        /*
        orsa::Matrix m;
        if (!orsa::eulerAnglesToMatrix(m,
                                       par->get("v2e_psi"),
                                       par->get("v2e_theta"),
                                       par->get("v2e_phi"))) {
            ORSA_DEBUG("problems...");
            exit(0);
        }
        const orsa::Matrix v2e = m;
        if (!orsa::Matrix::invert(v2e,m)) {
            ORSA_DEBUG("problems...");
            exit(0);
        }
        const orsa::Matrix e2v = m;
        const double z_axis_angle = acos(orsa::Vector(0,0,1)*e2v*orsa::Vector(0,0,1));
        const orsa::Vector z_v = e2v*orsa::Vector(0,0,1);
        const double z_v_lat = asin(z_v.getZ());
        const double z_v_lon = atan2(z_v.getY(),z_v.getX());
        */
        /*********************************************/
        
        /*
        ORSA_DEBUG("a,b,c: %7.3f %7.3f %7.3f [km]   phi,theta,psi: %+8.3f %+8.3f %+8.3f [deg]   z-axis angle,lat,lon: %+8.3f %+8.3f %+8.3f [deg]   dv: %+7.3f %+7.3f %+7.3f [km]   RMS: %g [km]",
                   orsa::FromUnits(par->get("a"),orsa::Unit::KM,-1),
                   orsa::FromUnits(par->get("b"),orsa::Unit::KM,-1),
                   orsa::FromUnits(par->get("c"),orsa::Unit::KM,-1),
                   par->get("v2e_phi")*orsa::radToDeg(),
                   par->get("v2e_theta")*orsa::radToDeg(),
                   par->get("v2e_psi")*orsa::radToDeg(),
                   z_axis_angle*orsa::radToDeg(),
                   z_v_lat*orsa::radToDeg(),
                   z_v_lon*orsa::radToDeg(),
                   orsa::FromUnits(par->get("v2e_dx"),orsa::Unit::KM,-1),
                   orsa::FromUnits(par->get("v2e_dy"),orsa::Unit::KM,-1),
                   orsa::FromUnits(par->get("v2e_dz"),orsa::Unit::KM,-1),
                   orsa::FromUnits(retVal,orsa::Unit::KM,-1));
        */
        ORSA_DEBUG("a,b,c: %7.3f %7.3f %7.3f [km]   RMS: %g [km]",
                   orsa::FromUnits(par->get("a"),orsa::Unit::KM,-1),
                   orsa::FromUnits(par->get("b"),orsa::Unit::KM,-1),
                   orsa::FromUnits(par->get("c"),orsa::Unit::KM,-1),
                   orsa::FromUnits(retVal,orsa::Unit::KM,-1));
                   
        return retVal;
    }
protected:
    /* 
       void singleIterationDone(const orsa::MultiminParameters * par) const {
       ORSA_DEBUG("a: %8.3f   b: %8.3f   c: %8.3f",
       orsa::FromUnits(par->get("a"),orsa::Unit::KM,-1),
       orsa::FromUnits(par->get("b"),orsa::Unit::KM,-1),
       orsa::FromUnits(par->get("c"),orsa::Unit::KM,-1));
       }
    */
};

int main (int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if (argc != 4) {
        ORSA_DEBUG("Usage: %s <geoid.xyz> <lat-min-deg> <lat-max-deg>",argv[0]);
        exit(0);
    }
    
    FILE * fp = fopen(argv[1],"r");
    if (fp == 0) {
        ORSA_DEBUG("cannot open file [%s]",argv[1]);
        exit(0);
    }
    
    const double lat_min_deg = std::min(atof(argv[2]),atof(argv[3]));
    const double lat_max_deg = std::max(atof(argv[2]),atof(argv[3]));
    
    /* const double lat_min = orsa::degToRad()*lat_min_deg;
       const double lat_max = orsa::degToRad()*lat_max_deg;
    */
    
    GeoidFitMultimin::VertexVector vertexVector;
    
    double lon_deg,lat_deg,val;
    char line[1024];
    while (fgets(line,1024,fp) != 0) {
        if (3 == sscanf(line,"%lf %lf %lf",&lon_deg,&lat_deg,&val)) {
            vertex v;
            v.lon = lon_deg*orsa::degToRad();
            v.lat = lat_deg*orsa::degToRad();
            double s_lon, c_lon;
           	// sincos(v.lon,&s_lon,&c_lon);
            s_lon = sin(v.lon);
			c_lon = cos(v.lon);
			double s_lat, c_lat;
 			// sincos(v.lat,&s_lat,&c_lat);
            s_lat = sin(v.lat);
			c_lat = cos(v.lat);
			v.v = orsa::FromUnits(val,orsa::Unit::KM)*orsa::Vector(c_lon*c_lat,s_lon*c_lat,s_lat);
            if ( (lat_deg >= lat_min_deg) &&
                 (lat_deg <= lat_max_deg) ) {
                v.use = true;
            } else {
                v.use = false;
            }
            vertexVector.push_back(v);
        }
    }
    
    fclose (fp);
    
    osg::ref_ptr<GeoidFitMultimin> multimin =
        new GeoidFitMultimin(vertexVector);
    
    osg::ref_ptr<orsa::MultiminParameters> par =
        new orsa::MultiminParameters;
    par->insert("a",orsa::FromUnits(280.0,orsa::Unit::KM),orsa::FromUnits(10.0,orsa::Unit::KM));
    par->insert("b",orsa::FromUnits(275.0,orsa::Unit::KM),orsa::FromUnits(10.0,orsa::Unit::KM));
    par->insert("c",orsa::FromUnits(230.0,orsa::Unit::KM),orsa::FromUnits(10.0,orsa::Unit::KM));
    // v2e =  vertex ref. sys.  TO  ellipse ref. sys.
    /*
    par->insert("v2e_psi",   0.0*orsa::degToRad(), 90.0*orsa::degToRad());
    par->insert("v2e_theta", 0.0*orsa::degToRad(), 30.0*orsa::degToRad());
    par->insert("v2e_phi",   0.0*orsa::degToRad(), 90.0*orsa::degToRad());
    */
    //
    // center offset
    /*
    par->insert("v2e_dx",orsa::FromUnits(0.0,orsa::Unit::KM),orsa::FromUnits(3.0,orsa::Unit::KM));
    par->insert("v2e_dy",orsa::FromUnits(0.0,orsa::Unit::KM),orsa::FromUnits(3.0,orsa::Unit::KM));
    par->insert("v2e_dz",orsa::FromUnits(0.0,orsa::Unit::KM),orsa::FromUnits(3.0,orsa::Unit::KM));
    */
    
    multimin->setMultiminParameters(par.get());
    const bool converged = multimin->run_nmsimplex();
    // const bool converged = multimin->run_nmsimplex2();
    // const bool converged = multimin->run_conjugate_fr();
    
    if (!converged) {
        ORSA_DEBUG("problems...");
    } else {
        GeoidFitMultimin::VertexVector ellipsoidVertexVector;
        GeoidFitMultimin::getEllipsoidVV(ellipsoidVertexVector,vertexVector,par,true);
        FILE * fp = fopen("geoid_ref_ellipsoid.xyz","w");
        ORSA_DEBUG("writing file [geoid_ref_ellipsoid.xyz]...");
        for (size_t k=0; k<ellipsoidVertexVector.size(); ++k) {
            const double sign = (vertexVector[k].v-ellipsoidVertexVector[k].v)*ellipsoidVertexVector[k].v > 0.0 ? +1 : -1;
            gmp_fprintf(fp,"%g %g %g\n",
                       orsa::radToDeg()*vertexVector[k].lon,
                       orsa::radToDeg()*vertexVector[k].lat,
                       orsa::FromUnits(sign*(vertexVector[k].v-ellipsoidVertexVector[k].v).length(),orsa::Unit::KM,-1));
#warning plain distance or projected distance?
            
#warning THIS should change: output plain ellipsoid in xyz, and then use GMT to subtract geoid and such...
            
        }
        fclose(fp);
    }
    
    return 0;
}
