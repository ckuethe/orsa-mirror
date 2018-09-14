#include <orsa/debug.h>
#include <cstdlib>
#include <cstdio>

#include <orsa/debug.h>
#include <orsa/statistic.h>
#include <orsa/util.h>
#include <orsa/vector.h>
#include <osg/ref_ptr>

class VecVal {
public:
    orsa::Vector u;
    double val;
};

void read(std::vector<VecVal> & VV, const std::string & filename, const double & c_over_a) {
    FILE * fp = fopen(filename.c_str(),"r");
    double lon_deg,lat_deg,val;
    char line[1024];
    while (fgets(line,1024,fp) != 0) {
        if (3 == sscanf(line,"%lf %lf %lf",&lon_deg,&lat_deg,&val)) {
            double lon = lon_deg*orsa::degToRad();
            double lat = lat_deg*orsa::degToRad();
            double s_lon, c_lon;
            s_lon = sin(lon);
            c_lon = cos(lon);
            double s_lat, c_lat;
            s_lat = sin(lat);
            c_lat = cos(lat);
            //
            VecVal vv;
            vv.u = orsa::Vector(c_lon*c_lat,s_lon*c_lat,s_lat*c_over_a).normalized();
            vv.val = val;
            VV.push_back(vv);
        }
    }
    fclose(fp);
}

double val(const std::vector<VecVal> & VV, const orsa::Vector & u) {
    double best_prod=-1.0;
    double best_val;
    double prod;
    for (size_t k=0; k<VV.size(); ++k) {
        prod=VV[k].u*u;
        if (prod>best_prod) {
            best_prod=prod;
            best_val=VV[k].val;
        }
    }
    return best_val;
}

int main (int argc, char **argv) {

    if (argc != 3) {
        ORSA_DEBUG("Usage: %s <file.xyz> <lon-lat-file>",argv[0]);
        exit(0);
    }
    
    const double c_over_a = 1.0;
    
    std::vector<VecVal> XYZ;
    read(XYZ,argv[1],c_over_a);
    
    FILE * fp = fopen(argv[2],"r");
    if (fp == 0) {
        ORSA_DEBUG("cannot open file [%s]",argv[2]);
        exit(0);
    }
    
    double lon_deg,lat_deg;
    char line[1024];
    while (fgets(line,1024,fp) != 0) {
        if (2 == sscanf(line,"%lf %lf",&lon_deg,&lat_deg)) {
            const double lon = lon_deg*orsa::degToRad();
            const double lat = lat_deg*orsa::degToRad();
            //
            double s_lon, c_lon;
            s_lon = sin(lon);
            c_lon = cos(lon);
            //
            double s_lat, c_lat;
            s_lat = sin(lat);
            c_lat = cos(lat);
            //
            const orsa::Vector u = orsa::Vector(c_lon*c_lat,s_lon*c_lat,s_lat*c_over_a).normalized();
            
            const double value = val(XYZ,u);
            printf("%g\n",value);
        }
    }
    
    fclose (fp);
        
    return 0;
}
