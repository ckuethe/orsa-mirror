#include "ThermalStressCraters.h"

int main (int argc, char **argv) {
    
    if (argc != 9999) {
        ORSA_DEBUG("Usage: %s <D,km> <d,km> <lat,deg> <lon,deg> <R,km> <slope-out,deg> <slope-azimuth,deg> <slope-rim,deg> <shape-par> <pole-R.A.,deg> <pole-Dec.,deg> <pole-phi-J2000,deg> <a,AU> <ecc> <i,deg> <node,deg> <peri,deg> <off-North,km> <off-East,km>",argv[0]);
        exit(0);
    }
    
    
    
    return 0;
}
