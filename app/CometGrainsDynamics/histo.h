#ifndef CGD_HISTO_H
#define CGD_HISTO_H

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

#endif // CGD_HISTO_H
