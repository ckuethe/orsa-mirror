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

class Mask {
public:
    orsa::Vector u; // center direction, unit vector
    // double phi; // angular distance from u
    double cos_phi; // cosine of angle of RADIUS of crater, not DIAMETER
};

bool test_masked(const orsa::Vector & u, const std::vector<Mask> & mask) {
    bool masked=false;
    for(size_t q=0; q<mask.size(); ++q) {
        if (u*mask[q].u > mask[q].cos_phi) {
            masked=true;
            break;
        }
    }
    return masked;
}

void read(std::vector<Mask> & mask, const std::string & filename) {
    FILE * fp = fopen(filename.c_str(),"r");
    double lon_deg,lat_deg,phi_deg;
    char line[1024];
    while (fgets(line,1024,fp) != 0) {
        if (line[0]=='#') continue;
        if (3 == sscanf(line,"%lf %lf %lf",&lon_deg,&lat_deg,&phi_deg)) {
            double lon = lon_deg*orsa::degToRad();
            double lat = lat_deg*orsa::degToRad();
            double phi = phi_deg*orsa::degToRad();
            double s_lon, c_lon;
            s_lon = sin(lon);
            c_lon = cos(lon);
            double s_lat, c_lat;
            s_lat = sin(lat);
            c_lat = cos(lat);
            double cos_phi = cos(phi);
            //
            Mask m;
            m.u = orsa::Vector(c_lon*c_lat,s_lon*c_lat,s_lat).normalized();
            // m.phi = val;
            m.cos_phi = cos_phi;
            mask.push_back(m);
        }
    }
    fclose(fp);
}

orsa::Vector sample_u(const double & c_over_a) {
    double dx,dy,dz;
    orsa::GlobalRNG::instance()->rng()->gsl_ran_dir_3d(&dx,&dy,&dz);
    const orsa::Vector u = orsa::Vector(dx,dy,dz*c_over_a).normalized();
    return u;    
}

orsa::Vector sample_u(const std::vector<Mask> & mask, const double & c_over_a) {
    orsa::Vector u;
    bool masked=true;
    while (masked) {
        u = sample_u(c_over_a);
        masked=test_masked(u,mask);
    }
    return u;
}

int main (int argc, char **argv) {

    /*
    if ((argc!=4) && (argc!=5) && (argc!=6)) {
        ORSA_DEBUG("Usage: %s <file1.xyz> <file2.xyz> <num-samples-per-iteration> [(c/a)=1.0] [mask.xyz]",argv[0]);
        exit(0);
    }
    */
    //
    /*
    if ((argc!=6) && (argc!=7) && (argc!=8)) {
        ORSA_DEBUG("Usage: %s <X.xyz> <X_mask.xyz> <Y.xyz> <Y_mask.xyz> <num-samples-per-iteration> [(c/a)=1.0] [random-rotation=0]",argv[0]);
        exit(0);
    }
    */
    //
    if ((argc!=7) && (argc!=8) && (argc!=9)) {
        ORSA_DEBUG("Usage: %s <X.xyz> <X_mask.xyz> <Y.xyz> <Y_mask.xyz> <num-samples-per-iteration> <max-iterations> [(c/a)=1.0] [random-rotation=0]",argv[0]);
        exit(0);
    }
    
    double c_over_a=1.0;
    if (argc>=8) {
        c_over_a = atof(argv[7]);
    }
    
    std::vector<VecVal> X;
    read(X,argv[1],c_over_a);
    
    std::vector<Mask> X_mask;
    read(X_mask,argv[2]);
    
    std::vector<VecVal> Y;
    read(Y,argv[3],c_over_a);
    
    std::vector<Mask> Y_mask;
    read(Y_mask,argv[4]);
    
    const size_t num_samples_per_iteration = atoi(argv[5]);
    const size_t max_iterations = atoi(argv[6]);
    
    // c/a above for argv[7]
    
    bool random_rotation=false;
    if (argc>=9) {
        random_rotation = atoi(argv[8]);
    }
    
    ORSA_DEBUG("X.size(): %u   X_mask.size(): %u   Y.size(): %u   Y_mask.size(): %u   samples: %i   iterations: %i   c/a: %g   rot: %i",X.size(),X_mask.size(),Y.size(),Y_mask.size(),num_samples_per_iteration,max_iterations,c_over_a,random_rotation);
    
    osg::ref_ptr<orsa::Statistic<double> > stat_rho = new orsa::Statistic<double>;
    
    const std::vector<Mask> empty_mask;
    
    FILE * fp_dump = fopen("correlation_dump.dat","w");
    
    size_t iter = 0;
    while (iter<max_iterations) {
        ++iter;
        
        osg::ref_ptr<orsa::Statistic<double> > stat_X  = new orsa::Statistic<double>;
        osg::ref_ptr<orsa::Statistic<double> > stat_XX = new orsa::Statistic<double>;
        osg::ref_ptr<orsa::Statistic<double> > stat_Y  = new orsa::Statistic<double>;
        osg::ref_ptr<orsa::Statistic<double> > stat_YY = new orsa::Statistic<double>;
        osg::ref_ptr<orsa::Statistic<double> > stat_XY = new orsa::Statistic<double>;
        
        // relative rotation of Y with respect to X to test significance by assigning random rotations
        orsa::Matrix rot = orsa::Matrix::identity();
        if (random_rotation) {
            orsa::Vector u = sample_u(empty_mask,c_over_a);
            const double phi   = atan2(u.getY(),u.getX());
            const double theta = acos(u.getZ());
            const double psi   = orsa::twopi()*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
            orsa::eulerAnglesToMatrix(rot,psi,theta,phi);
        }
        double last_rho;
        for (size_t jj=0; jj<num_samples_per_iteration; ++jj) {
            orsa::Vector u_X;
            orsa::Vector u_Y;
            bool good_uu=false;
            while (good_uu==false) {
                u_X = sample_u(c_over_a);
                u_Y = rot*u_X;
                good_uu=true;
                if (test_masked(u_X,X_mask)) {
                    good_uu=false;
                    continue;
                }
                if (test_masked(u_Y,Y_mask)) {
                    good_uu=false;
                    continue;
                }
            }
            const double X_val = val(X,u_X);
            const double Y_val = val(Y,u_Y);
            {
                // debug only...
                const double X_lat = asin(u_X.getZ());
                const double X_lon = atan2(u_X.getY(),u_X.getX());
                const double Y_lat = asin(u_Y.getZ());
                const double Y_lon = atan2(u_Y.getY(),u_Y.getX());
                fprintf(fp_dump,"%g %g %g %g %g %g\n",fmod(orsa::twopi()+X_lon,orsa::twopi())*orsa::radToDeg(),X_lat*orsa::radToDeg(),fmod(orsa::twopi()+Y_lon,orsa::twopi())*orsa::radToDeg(),Y_lat*orsa::radToDeg(),X_val,Y_val);
            }
            stat_X->insert(X_val);
            stat_Y->insert(Y_val);
            stat_XX->insert(X_val*X_val);
            stat_YY->insert(Y_val*Y_val);
            stat_XY->insert(X_val*Y_val);
        }
        
        // expectation values
        const double E_X  = stat_X->average();
        const double E_Y  = stat_Y->average();
        const double E_XX = stat_XX->average();
        const double E_YY = stat_YY->average();
        const double E_XY = stat_XY->average();
        
        // correlation coefficient rho
        const double rho = (E_XY-E_X*E_Y)/(sqrt(E_XX-E_X*E_X)*sqrt(E_YY-E_Y*E_Y));
        last_rho=rho;
        stat_rho->insert(last_rho);
        
        ORSA_DEBUG("rho: %.3f +/- %.3f   stddev: %.3f",stat_rho->average(),stat_rho->averageError(),stat_rho->standardDeviation());
        
        FILE * fp_latest = fopen("correlation_latest.dat","w");
        fprintf(fp_latest,"%+.3f %.3f %.3f\n",stat_rho->average(),stat_rho->averageError(),stat_rho->standardDeviation());
        fclose(fp_latest);
    }
    
    fclose (fp_dump);
    
    return 0;
}
