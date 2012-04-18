#include "histo.h"

int main (int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    const double pixelScale = orsa::FromUnits(20.0,orsa::Unit::KM);
    // every axis has 2*N pixels, from -N to N-1
    const size_t N = 100;
    
    // const double min_dt = 
    
    FILE * fp = fopen("colden.out","r");
    if (!fp) {
        exit(0);
    }
    
    std::vector<ColDenData> colden;
    {
        char line[4096];
        double r0,rF,AF,dt,pos_sun,pos_orbit_pole,pos_orbit_plane;
#warning keep fields in sync with CometGrainsDynamics.cpp
        while (fgets(line,4096,fp)) {
            gmp_sscanf(line,
                       "%*s   %*s %*s   %lf %*s   %lf %lf %lf   %*s   %lf %lf %lf   %*s",
                       &dt,
                       &r0,
                       &rF,
                       &AF,
                       &pos_sun,
                       &pos_orbit_pole,
                       &pos_orbit_plane);
            dt = orsa::FromUnits(dt,orsa::Unit::DAY);
            r0 = orsa::FromUnits(r0,orsa::Unit::METER);
            rF = orsa::FromUnits(rF,orsa::Unit::METER);
            AF = orsa::FromUnits(AF,orsa::Unit::METER,2);
            pos_sun = orsa::FromUnits(pos_sun,orsa::Unit::KM);
            pos_orbit_pole = orsa::FromUnits(pos_orbit_pole,orsa::Unit::KM);
            pos_orbit_plane = orsa::FromUnits(pos_orbit_plane,orsa::Unit::KM);
            ColDenData data;
            data.r0 = r0;
            data.rF = rF;
            data.AF = AF;
            data.dt = dt;
            data.pos_sun = pos_sun;
            data.pos_orbit_pole = pos_orbit_pole;
            data.pos_orbit_plane = pos_orbit_plane;
            colden.push_back(data);
        }
        ORSA_DEBUG("entries: %i",colden.size());
        fclose(fp);
    }

    {
        // 1D 
        std::vector<double> histo_1D;
        histo_1D.resize(2*N);
        for (size_t k=0; k<2*N; ++k) {
            histo_1D[k] = 0.0;
        }
        
        for (size_t j=0; j<colden.size(); ++j) {
            const double  pos = colden[j].pos_sun;
            const double area = colden[j].AF;
            if (fabs(pos)>N*pixelScale) continue;
            histo_1D[N+pos/pixelScale] += area;
        }
        
        for (size_t k=0; k<2*N; ++k) {
            ORSA_DEBUG("%i %g %g",
                       k,
                       orsa::FromUnits((k+0.5-N)*pixelScale,orsa::Unit::KM,-1),
                       orsa::FromUnits(histo_1D[k],orsa::Unit::CM,-2));
        }
    }
    
    {
        // 2D
        
        
    }
    
    return 0;
}
