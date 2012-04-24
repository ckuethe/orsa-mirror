#include "histo.h"

#include <algorithm>
#include <deque>
#include <list>

int main (int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if (argc != 2) {
        ORSA_DEBUG("Usage: %s <input-colden-file>",argv[0]);
        exit(0);
    }
    
    const std::string inputFile = argv[1];
    
    const double pixelScale = orsa::FromUnits(50.0,orsa::Unit::KM);
    // every axis has 2*N pixels, from -N to N-1
    const int N = __ORSA_FITS_N__/2;
    
    // keep this in sync with CometGrainsDynamics.cpp
    const double min_dt = orsa::FromUnits(  5,orsa::Unit::SECOND); // cannot be zero!
    const double max_dt = orsa::FromUnits(200,orsa::Unit::DAY);
    //
    const double min_grainRadius = orsa::FromUnits(1.0e-6,orsa::Unit::METER);
    const double max_grainRadius = orsa::FromUnits(1.0   ,orsa::Unit::METER);    
    // population scaling law
    const double pop_N0 = 1.0;
    const double pop_a  = 4.5; // this is positive, the minus sign is already included in code
    
    const size_t min_selected = 1; // 2
    const size_t min_grains_per_pixel = 1; // 2
    
    FILE * fp = fopen(inputFile.c_str(),"r");
    if (!fp) {
        ORSA_DEBUG("cannot open file [%s]",inputFile.c_str());
        exit(0);
    }
    
    typedef std::list<ColDenData> ColDenDataContainer;
    ColDenDataContainer colden;
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
        // ORSA_DEBUG("entries: %i",colden.size());
        fclose(fp);
    }
    const size_t colden_size = colden.size();
    
    // important for the lower_bound(...) calls below
    colden.sort();
    //
    /* 
       {
       // test
       ColDenDataContainer::iterator it = colden.begin();
       while (it != colden.end()) {
       ORSA_DEBUG("dt: %g",(*(*it).dt));
       ++it;
       }
       }
    */
    
    // colden is destroyed by the algorithm following, so if reuse is required,
    // it must be copied to another variable first
    
    {        
        // 2D
        std::vector< std::vector< std::deque<double> > > histo_2D;
        histo_2D.resize(2*N);
        for (size_t i=0; i<2*N; ++i) {
            histo_2D[i].resize(2*N);
        }
        
        size_t zeroArea=0;
        size_t inField=0;
        
        // two outer loops: one on the time interval, one on the grain radius
        
        ColDenData colden_lower_bound_dt_value;
        
        double  low_dt = min_dt;
        double high_dt;
        while (low_dt < max_dt) {
            high_dt = std::min(1.50*low_dt,max_dt);
            // high_dt = std::min(low_dt+orsa::FromUnits(1.0,orsa::Unit::HOUR),max_dt);
            // ORSA_DEBUG("dt range: %g %g",low_dt,high_dt);
            
            colden_lower_bound_dt_value.dt = low_dt;
            ColDenDataContainer::iterator it_low_dt =
                lower_bound(colden.begin(), colden.end(), colden_lower_bound_dt_value);
            
            colden_lower_bound_dt_value.dt = high_dt;
            ColDenDataContainer::iterator it_high_dt =
                lower_bound(colden.begin(), colden.end(), colden_lower_bound_dt_value);
            
            // ORSA_DEBUG("dt: %g it: %g", low_dt, (*(*it_low_dt).dt));
            // ORSA_DEBUG("dt: %g it: %g",high_dt,(*(*it_high_dt).dt));
            
            if (it_low_dt == it_high_dt) {
                low_dt = high_dt;
                continue;
            }
            
            ColDenDataContainer selected_dt_colden;
            selected_dt_colden.splice(selected_dt_colden.begin(),
                                      colden,
                                      it_low_dt,
                                      it_high_dt);
            /*
               {
               // test
               ColDenDataContainer::iterator it = selected_dt_colden.begin();
               while (it != selected_dt_colden.end()) {
               ORSA_DEBUG("dt: %g",(*(*it).dt));
               ++it;
               }
               }
            */
            
            double  low_Rg = min_grainRadius;
            double high_Rg = low_Rg;
            bool   used_Rg = false;
            while (1) {
                high_Rg = std::min(1.20*high_Rg,max_grainRadius); // self-referring
                // ORSA_DEBUG("Rg range: %g %g",low_Rg,high_Rg);
                ColDenDataContainer selected_dt_Rg_colden;
                ColDenDataContainer::iterator it = selected_dt_colden.begin();
                while (it != selected_dt_colden.end()) {
                    /* if ( ((*it).dt >= low_dt) &&
                       ((*it).dt < high_dt) &&
                       ((*it).r0 >= low_Rg) &&
                       ((*it).r0 < high_Rg) ) {
                    */
                    if ( // ((*it).dt >= low_dt) &&
                         // ((*it).dt < high_dt) &&
                        ((*it).r0 >= low_Rg) &&
                        ((*it).r0 < high_Rg) ) {
                        selected_dt_Rg_colden.push_back((*it));
                    }
                    ++it;
                }
                // ORSA_DEBUG("selected for this bin: %i",selected_colden.size());
                // if (selected_colden.size()>0) {
                
                /* ORSA_DEBUG("dt: [%.3e-%.3e]   Rg: [%.3e-%.3e]   pop: %.3e   selected: %i",
                   low_dt,high_dt,
                   low_Rg,high_Rg,
                   pop(low_Rg,high_Rg,pop_N0,pop_a),
                   selected_colden.size());
                */
                
                if (selected_dt_Rg_colden.size()==0) {
                    used_Rg=true;
                } else if (selected_dt_Rg_colden.size()>=min_selected) {                 
                    const double factor = (high_dt-low_dt)*pop(low_Rg,high_Rg,pop_N0,pop_a)/selected_dt_Rg_colden.size();
                    /* ORSA_DEBUG("dt: [%.3e-%.3e]   Rg: [%.3e-%.3e]   factor: %.3e   pop: %.3e   selected: %i",
                       low_dt,high_dt,
                       low_Rg,high_Rg,
                       factor,
                       pop(low_Rg,high_Rg,pop_N0,pop_a),
                       selected_colden.size());
                    */
                    // for (size_t l=0; l<selected_colden.size(); ++l) {
                    ColDenDataContainer::iterator it = selected_dt_Rg_colden.begin();
                    while (it != selected_dt_Rg_colden.end()) {
                        if ((*it).AF==0.0) {
                            ++zeroArea;
                            ++it;
                            continue;
                        }
                        //
                        const double  pos_i = (*it).pos_sun;
                        const int i = N+pos_i/pixelScale;
                        if (i<0) {
                            ++it;
                            continue;
                        }
                        if (i>=2*N) {
                            ++it;
                            continue; 
                        }
                        //
                        const double  pos_j = (*it).pos_orbit_plane;
                        const int j = N+pos_j/pixelScale;
                        if (j<0) {
                            ++it;
                            continue;
                        }
                        if (j>=2*N) {
                            ++it;
                            continue;
                        }
                        //
                        const double scaled_area = factor*(*it).AF;
                        
                        // choose one
                        //
                        histo_2D[i][j].push_back(scaled_area);
                        //
                        // rho (projected distance from nucleus) is for testing 1/rho
                        /* const double rho = sqrt(orsa::square(pos_i)+orsa::square(pos_j));
                           histo_2D[i][j].push_back(scaled_area*rho);
                        */
                        
                        ++inField;
                        ++it;
                    }
                    
                    used_Rg=true;
                }
                
                if (used_Rg) {
                    low_Rg = high_Rg;
                    used_Rg=false;
                }
                
                if (high_Rg>=max_grainRadius) break;
            }
            
            low_dt = high_dt;
        }

        size_t populatedPixels=0;
        
        // #warning write zeroes?
        FILE * fp = fopen("histo_2D.out","w");
        for (size_t i=0; i<2*N; ++i) {
            for (size_t j=0; j<2*N; ++j) {
                if (histo_2D[i][j].size()>=min_grains_per_pixel) {
                    double val = 0.0;
                    for (size_t k=0; k<histo_2D[i][j].size(); ++k) {
                        val += histo_2D[i][j][k];
                    }
                    if (val>0.0) {
                        ++populatedPixels;
                        fprintf(fp,"%3i %3i %.3e\n",i,j,val);
                    }
                }
            }
        }
        fclose(fp);
        
        ORSA_DEBUG("entries: %i   in-field: %i   zero-area: %i   pixels: %i (%.1f\%)",colden_size,inField,zeroArea,populatedPixels,100.0*(double)populatedPixels/(4*N*N));
        
    }
    
    return 0;
}
