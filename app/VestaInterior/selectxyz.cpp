#include <orsa/debug.h>
#include <cstdlib>
#include <cstdio>

int main (int argc, char **argv) {

    if (argc != 4) {
        ORSA_DEBUG("Usage: %s <file.xyz> <lat-min-deg> <lat-max-deg>",argv[0]);
        exit(0);
    }
    
    FILE * fp = fopen(argv[1],"r");
    if (fp == 0) {
        ORSA_DEBUG("cannot open file [%s]",argv[1]);
        exit(0);
    }
    
    const double lat_min_deg = std::min(atof(argv[2]),atof(argv[3]));
    const double lat_max_deg = std::max(atof(argv[2]),atof(argv[3]));

    double lon_deg,lat_deg,val;
    char line[1024];
    while (fgets(line,1024,fp) != 0) {
        if (3 == sscanf(line,"%lf %lf %lf",&lon_deg,&lat_deg,&val)) {
            if ( (lat_deg >= lat_min_deg) &&
                 (lat_deg <= lat_max_deg) ) {
                printf("%s",line);
            }
        }   
    }
    
    fclose (fp);
    
    return 0;
}
