#include "ThermalStressCraters.h"

int main (int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if (argc != 9999) {
        ORSA_DEBUG("Usage: %s <D,km> <d,km> <lat,deg> <lon,deg> <R,km> <slope-out,deg> <slope-azimuth,deg> <slope-rim,deg> <shape-par> <pole-R.A.,deg> <pole-Dec.,deg> <pole-phi-J2000,deg> <a,AU> <ecc> <i,deg> <node,deg> <peri,deg> <off-North,km> <off-East,km>",argv[0]);
        // exit(0);
    }


    if (0) {
        // test
        double r=0.0;
        double h,dhdr;
        while (r<=50.0) {
            CraterShape(h,dhdr,r,100,10,tan(0.0*orsa::degToRad()),tan(40.0*orsa::degToRad()));
            ORSA_DEBUG("cs %g %g %g",r,h,atan(dhdr)*orsa::radToDeg());
            r+=0.01;
        }
    }
    
    const size_t numSlices=400;
    History history;
    const size_t NS=100000;
    const size_t days=70;
    const double hdist = orsa::FromUnits(3.0,orsa::Unit::AU); // heliocentric distance
    ORSA_DEBUG("big-theta: %g",theta(hdist));
    std::vector<double> Fs;
    Fs.resize(NS);
    for (size_t p=0; p<NS; ++p) {
        // #warning restore this one!
        // Fs[p] = solar()/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2)*std::max(cos(days*orsa::twopi()*(double)p/(double)NS),0.0);

        // test
        Fs[p] = solar()/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2) *
            std::max(cos(days*orsa::twopi()*(double)p/(double)NS),0.0) *
            std::max(0.5+0.7*cos(orsa::twopi()*(double)p/(double)NS),0.0);
        if (p%99==0) Fs[p] = 0.0;
        
        // TEST!
        /* double proj = std::max(cos(days*orsa::twopi()*(double)p/(double)NS),0.0);
           if (proj < 0.5) proj = 0.0;
           Fs[p] = solar()/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2)*proj;
        */
    }
    
    /* for (size_t p=0; p<NS; ++p) {
       ORSA_DEBUG("Fs[%06i] = %g",p,Fs[p]);
       }
    */
    
    const double dx = 0.2*skinDepth(); // or a fraction of skinDepth = ls
    const double dt = days*rotationPeriod()/NS;
    const unsigned int history_skip = 1;
    
    ComputePeriodicThermalHistory(history,
                                  numSlices,
                                  200.0,
                                  Fs,
                                  dx,
                                  dt,
                                  history_skip);
    
    
    FILE * fp = fopen("TSC.out","w");
    for (unsigned int j=0; j<history.size(); ++j) {

        /* unsigned int jmm = (j==0) ? (history.size()-1) : (j-1);
           const double dTdt = (history[j][0].T-history[jmm][0].T)/dt;
        */
        
        unsigned int jpp = (j==(history.size()-1)) ? (0) : (j+1);
        const double dTdt = (history[jpp][0].T-history[j][0].T)/dt;
        
        gmp_fprintf(fp,
                    "%g %g %g %g %g\n",
                    j*history_skip*dt,
                    Fs[j*history_skip],
                    history[j][0].T,
                    history[j][history[j].size()-2].T, // -2 because -1 is never changed by thermal algo
                    dTdt);
        
    }
    fclose(fp);
    
    return 0;
}
