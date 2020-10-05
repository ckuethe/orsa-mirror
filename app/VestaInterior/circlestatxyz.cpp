#include <orsa/debug.h>
#include <orsa/multimin.h>
#include <orsa/print.h>
#include <orsa/shape.h>
#include <orsa/statistic.h>
#include <orsa/util.h>
#include <orsa/vector.h>

class VecVal {
public:
    orsa::Vector u;
    double val;
};

/*
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
*/

void read(std::vector<VecVal> & VV, const std::string & filename) {
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
            vv.u = orsa::Vector(c_lon*c_lat,s_lon*c_lat,s_lat).normalized();
            vv.val = val;
            VV.push_back(vv);
        }
    }
    fclose(fp);
}

/*
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
*/

int main (int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
        
    /*
    if (argc != 3) {
        ORSA_DEBUG("Usage: %s <file.xyz> <band-half-width-DEG>",argv[0]);
        exit(0);
    }
    */
    //
    if (argc != 4) {
        ORSA_DEBUG("Usage: %s <file.xyz> <band-half-width-DEG> <num-sample-points>",argv[0]);
        exit(0);
    }
    
    // c_over_a not needed since working in lat/lon in both input and output... right?
    /*
    const double c_over_a = 0.9253;
    ORSA_DEBUG("c_over_a: %g",c_over_a);
    */
    
    std::vector<VecVal> XYZ;
    // read(XYZ,argv[1],c_over_a);
    read(XYZ,argv[1]);
    const double band_half_width = atof(argv[2])*orsa::degToRad();
    const int max_num = atoi(argv[3]);
    
    int count=0;
    while (count<max_num) {
        
        // pole
        // should take into account polar flattenig of body...
        /*
        const double lon0 = orsa::twopi()*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
        const double tilt = 90.0*orsa::degToRad()*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
        */
        // use this to generate uniformly
        double x,y,z;
        orsa::GlobalRNG::instance()->rng()->gsl_ran_dir_3d(&x,&y,&z);
        if (z<0.0) {
            // skip: one hemisphere is enough...
            continue;
        }
        //
        const orsa::Vector u_pole = orsa::Vector(x,y,z).normalized();
        // lon0 is the longitude where the test-equator crosses the current equator, going up
        const double lon0 = fmod(atan2(y,x)+orsa::halfpi()+orsa::twopi(),orsa::twopi());
        const double tilt = orsa::halfpi()-fabs(asin(z));
        
        // NEW version: points within delta_DEG from test-equator
        osg::ref_ptr<orsa::Statistic<double> > stat   = new orsa::Statistic<double>;
        osg::ref_ptr<orsa::Statistic<double> > stat_w = new orsa::Statistic<double>;
        //
        const double max_SP = cos(orsa::halfpi()-band_half_width); // maximum value of scalar product of points with test-pole, 
        const double min_SP = cos(orsa::halfpi()+band_half_width); // minimum value of scalar product of points with test-pole, 
        //
        std::vector<VecVal>::const_iterator it = XYZ.begin();
        while (it != XYZ.end()) {
            const double SP = u_pole*(*it).u;
            // ORSA_DEBUG("%g %g    %g",min_SP,max_SP,SP);
            if ((SP>min_SP) && (SP<max_SP)) {
                const double w = 1.0; // isotropic grid
                // const double w = (1.0-orsa::square((*it).u.getZ())); // use this for regular lat/lon grid data
                stat->insert(w*(*it).val);
                stat_w->insert(w);
                // ORSA_DEBUG("%g %g %g",(*it).u.getX(),(*it).u.getY(),(*it).u.getZ());
            }
            ++it;
        }
        // printf("%g %g %g %g   %lu\n",lon0*orsa::radToDeg(),tilt*orsa::radToDeg(),stat->average(),stat->averageError(),stat->entries().get_ui());
        printf("%g %g %g %g\n",lon0*orsa::radToDeg(),tilt*orsa::radToDeg(),stat->sum()/stat_w->sum(),stat_w->sum());
        fflush(stdout);
        ++count;
    }
    
    return 0;
}
