#include <orsa/debug.h>
#include <orsa/double.h>
#include <orsa/util.h>
#include <cstdlib>
#include <cstdio>

class P {
public:
    double x,y,z;
};

int main (int argc, char **argv) {
    
    if (argc != 8) {
        ORSA_DEBUG("Usage: %s <file.xyz> <xmin> <xmax> <ymin> <ymax> <zmin> <zmax>",argv[0]);
        exit(0);
    }
    
    const double xmin = atof(argv[2]);
    const double xmax = atof(argv[3]);
    const double ymin = atof(argv[4]);
    const double ymax = atof(argv[5]);
    const double zmin = atof(argv[6]);
    const double zmax = atof(argv[7]);
    
    FILE * fp = fopen(argv[1],"r");
    if (fp == 0) {
        ORSA_DEBUG("cannot open file [%s]",argv[1]);
        exit(0);
    }
    
    std::vector<P> vec;
    vec.reserve(1000000);
    char line[1024];
    while (fgets(line,1024,fp) != 0) {
        P p;
        if (3 == sscanf(line,"%lf %lf %lf",&p.x,&p.y,&p.z)) {
            vec.push_back(p);
        }   
    }
    fclose (fp);

    ORSA_DEBUG("vec.size(): %i",vec.size());

    const double unset = -1.0;
    
    /* double d2xy_max=0.0;
       for (size_t k=0; k<vec.size(); ++k) {
       const P & pk = vec[k];
       double d2xy_k_min=unset;
       for (size_t j=0; j<k; ++j) {
       const P & pj = vec[j];
       const double d2xy = orsa::square(pj.x-pk.x)+orsa::square(pj.y-pk.y);
       if ( (d2xy < d2xy_k_min) ||
       (d2xy_k_min == unset) ) d2xy_k_min = d2xy;
       }
       if (d2xy_k_min > d2xy_max) {
       d2xy_max = d2xy_k_min;
       ORSA_DEBUG("new d2xy_max: %g",d2xy_max);
       }
       }
    */
    
    /* double dx_max=0.0;
       for (size_t k=0; k<vec.size(); ++k) {
       const P & pk = vec[k];
       if (pk.x < xmin) continue;
       if (pk.x > xmax) continue;
       double dx_k_min=unset;
       P p_k_min;
       for (size_t j=0; j<k; ++j) {
       const P & pj = vec[j];
       if (pj.x < xmin) continue;
       if (pj.x > xmax) continue;
       const double dx = fabs(pj.x-pk.x);
       if ( (dx < dx_k_min) ||
       (dx_k_min == unset) ) {
       dx_k_min = dx;
       p_k_min = pk;
       }
       }
       if (dx_k_min > dx_max) {
       dx_max = dx_k_min;
       ORSA_DEBUG("dx_max: %g   P: %g %g",dx_max,p_k_min.x,p_k_min.y);
       }
       }
    */

    /* 
       double dy_max=0.0;
       for (size_t k=0; k<vec.size(); ++k) {
       const P & pk = vec[k];
       if (pk.y < ymin) continue;
       if (pk.y > ymax) continue;
       double dy_k_min=unset;
       P p_k_min;
       for (size_t j=0; j<k; ++j) {
       const P & pj = vec[j];
       if (pj.y < ymin) continue;
       if (pj.y > ymax) continue;
       const double dy = fabs(pj.y-pk.y);
       if ( (dy < dy_k_min) ||
       (dy_k_min == unset) ) {
       dy_k_min = dy;
       p_k_min = pk;
       }
       }
       if (dy_k_min > dy_max) {
       dy_max = dy_k_min;
       ORSA_DEBUG("dy_max: %g   P: %g %g",dy_max,p_k_min.x,p_k_min.y);
       }
       }
    */
    
    /* double dx_max=0.0;
       double dy_max=0.0;
       for (size_t k=0; k<vec.size(); ++k) {
       const P & pk = vec[k];
       if (pk.x < xmin) continue;
       if (pk.x > xmax) continue;
       if (pk.y < ymin) continue;
       if (pk.y > ymax) continue;
       double dx_k_min=unset;
       double dy_k_min=unset;
       P p_k_min;
       for (size_t j=0; j<k; ++j) {
       const P & pj = vec[j];
       if (pj.x < xmin) continue;
       if (pj.x > xmax) continue;
       if (pj.y < ymin) continue;
       if (pj.y > ymax) continue;
       const double dx = fabs(pj.x-pk.x);
       const double dy = fabs(pj.y-pk.y);
       if ( ( (dx < dx_k_min) && (dy < dy_k_min) ) ||
       (dx_k_min == unset) ||
       (dy_k_min == unset) ) {
       dx_k_min = dx;
       dy_k_min = dy;
       p_k_min = pk;
       }
       }
       if ( (dx_k_min > dx_max) &&
       (dy_k_min > dy_max) ) {
       dx_max = dx_k_min;
       dy_max = dy_k_min;
       ORSA_DEBUG("dx_max: %g   dy_max: %g   P: %g %g",dx_max,dy_max,p_k_min.x,p_k_min.y);
       }
       }
    */
    
    /* while (1) {
       P p;
       p.x = xmin + (xmax-xmin)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
       p.y = ymin + (ymax-ymin)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
       bool good=true;
       for (size_t k=0; k<vec.size(); ++k) {
       const P & pk = vec[k];
       const double d2xy = orsa::square(p.x-pk.x)+orsa::square(p.y-pk.y);
       if (d2xy < d2xy_max) {
       good=false;
       break;
       }
       }
       if (good) {
       p.z = zmin + (zmax-zmin)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
       gmp_printf("%g %g %g\n",p.x,p.y,p.z);
       fflush(stdout);
       }
       }
    */

    // test: force values for dx_max and dy_max;
    const double dx_max = 0.003;
    const double dy_max = 0.20;
    
    while (1) {
        P p;
        p.x = xmin + (xmax-xmin)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
        p.y = ymin + (ymax-ymin)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
        bool good=true;
        for (size_t k=0; k<vec.size(); ++k) {
            const P & pk = vec[k];
            const double dx = fabs(p.x-pk.x);
            const double dy = fabs(p.y-pk.y);
            /* if ( (dx < dx_max) &&
               (dy < dy_max) ) {
               good=false;
               break;
               }
            */
            if ( orsa::square(dx/dx_max)+orsa::square(dy/dy_max) < 1.0 ) {
                good=false;
                break;
            }
        }
        if (good) {
            p.z = zmin + (zmax-zmin)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
            gmp_printf("%g %g %g\n",p.x,p.y,p.z);
            fflush(stdout);
        }
    }
    
    return 0;
}
