#ifndef CGD_HISTO_H
#define CGD_HISTO_H

#include "histosize.h"

#include <orsa/cache.h>
#include <orsa/datetime.h>
#include <orsa/debug.h>

class ColDenData {
public:
    orsa::Cache<double> r0; // initial grain radius (at t0)
    orsa::Cache<double> rF; // final grain radius (at t_snapshot)
    orsa::Cache<double> AF; // final area
    orsa::Cache<double> dt; // time between emission and snapshot
    orsa::Cache<double> pos_sun; 
    orsa::Cache<double> pos_orbit_pole;
    orsa::Cache<double> pos_orbit_plane;
};


// grain size distribution 
// differential distribution:
// N=N0*x^(-a), with x grain radius
// integrated between x1 and x2 gives:
// N2-N1 = N0/(a-1)*(x2^(-(a-1))-x1^(-(a-1)))
//
// also, this is per unit time
//
inline double pop(const double & minGrainRadius,
                  const double & maxGrainRadius,
                  const double & N0,
                  const double & a) {
    if (minGrainRadius>maxGrainRadius) {
        return pop(maxGrainRadius,minGrainRadius,N0,a);
    } else {
        return (N0/(1.0-a))*(pow(maxGrainRadius,1.0-a)-pow(minGrainRadius,1.0-a));
    }
}

#endif // CGD_HISTO_H
