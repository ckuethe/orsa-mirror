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

typedef std::vector<vertex> VertexVector;

void getEllipsoidVV(VertexVector & ellipsoidVertexVector, const VertexVector & vertexVector, const double & a, const double & b, const double & c, const bool useAll=false) {
    
    osg::ref_ptr<orsa::EllipsoidShape> shape = new orsa::EllipsoidShape(a,b,c);
    
    ellipsoidVertexVector.resize(vertexVector.size());
    orsa::Vector intersectionPoint;
    orsa::Vector normal;
    for (size_t k=0; k<vertexVector.size(); ++k) {
        ellipsoidVertexVector[k].use = vertexVector[k].use;
        if (useAll || vertexVector[k].use) {
            
            if (!shape->rayIntersection(intersectionPoint,normal,orsa::Vector(0,0,0),(*vertexVector[k].v).normalized(),false)) {
                ORSA_DEBUG("problems...");
                exit(0);
            }
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

int main (int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    const double km = orsa::FromUnits(1,orsa::Unit::KM);
    
    if (argc != 5) {
        ORSA_DEBUG("Usage: %s <file.xyz> <a/km> <b/km> <c/km>",argv[0]);
        exit(0);
    }
    
    FILE * fp = fopen(argv[1],"r");
    if (fp == 0) {
        ORSA_DEBUG("cannot open file [%s]",argv[1]);
        exit(0);
    }
    
    const double a = atof(argv[2])*km;
    const double b = atof(argv[3])*km;
    const double c = atof(argv[4])*km;
    
    VertexVector vertexVector;
    
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
            // if ( (lat_deg >= lat_min_deg) &&
            // (lat_deg <= lat_max_deg) ) {
                v.use = true;
            // } else {
            // v.use = false;
            // }
            vertexVector.push_back(v);
        }
    }
    fclose (fp);
    
    VertexVector ellipsoidVertexVector;
    getEllipsoidVV(ellipsoidVertexVector,vertexVector,a,b,c,true);
    FILE * fp_out = fopen("sub.xyz","w");
    ORSA_DEBUG("writing file [sub.xyz]...");
    for (size_t k=0; k<ellipsoidVertexVector.size(); ++k) {
        const double sign = (vertexVector[k].v-ellipsoidVertexVector[k].v)*ellipsoidVertexVector[k].v > 0.0 ? +1 : -1;
        gmp_fprintf(fp_out,"%g %g %g\n",
        orsa::radToDeg()*vertexVector[k].lon,
        orsa::radToDeg()*vertexVector[k].lat,
        orsa::FromUnits(sign*(vertexVector[k].v-ellipsoidVertexVector[k].v).length(),orsa::Unit::KM,-1));            
    }
    fclose(fp_out);

    return 0;
}
