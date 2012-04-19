#include "ThermalStressCraters.h"

int main (int argc, char **argv) {
    
    if (argc != 9999) {
        ORSA_DEBUG("Usage: %s <D,km> <d,km> <lat,deg> <lon,deg> <R,km> <slope-out,deg> <slope-azimuth,deg> <slope-rim,deg> <shape-par> <pole-R.A.,deg> <pole-Dec.,deg> <pole-phi-J2000,deg> <a,AU> <ecc> <i,deg> <node,deg> <peri,deg> <off-North,km> <off-East,km>",argv[0]);
        // exit(0);
    }



    const size_t NP=100;
    Profile profile(NP);
    const size_t NS=10000;
    const size_t days=7;
    std::vector<double> Fs;
    Fs.resize(NS);
    for (size_t p=0; p<NS; ++p) {
        Fs[p] = solar()*std::max(cos(days*orsa::twopi()*(double)p/(double)NS),0.0);
    }
    
    for (size_t p=0; p<NS; ++p) {
        ORSA_DEBUG("Fs[%06i] = %g",p,Fs[p]);
    }
    
    const double dx = skinDepth(); // or a fraction of skinDepth = ls
    const double dt = days*rotationPeriod()/NS;
    const unsigned int history_skip = 100;
    
    ComputePeriodicThermalProfile(profile,
                                  200.0,
                                  Fs,
                                  dx,
                                  dt,
                                  history_skip);
    
    return 0;
}
