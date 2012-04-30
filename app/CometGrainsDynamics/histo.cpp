#include "histo.h"

#include <algorithm>
#include <deque>
#include <list>

int main (int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if (argc == 1) {
        ORSA_DEBUG("Usage: %s <input-colden-file(s)>",argv[0]);
        exit(0);
    }
    
    const double pixelScale = orsa::FromUnits(1.0,orsa::Unit::KM);
    // every axis has 2*N pixels, from -N to N-1
    const int N = __ORSA_FITS_N__/2;
    
    // keep this in sync with CometGrainsDynamics.cpp
    const double min_dt = orsa::FromUnits(  5,orsa::Unit::SECOND); // cannot be zero!
    const double max_dt = orsa::FromUnits(200,orsa::Unit::DAY); // orsa::FromUnits(200,orsa::Unit::DAY);
    //
    const double min_grainRadius = orsa::FromUnits(0.010000,orsa::Unit::METER);
    const double max_grainRadius = orsa::FromUnits(1.000000,orsa::Unit::METER);    
    // population scaling law
    const double pop_N0 = 1.0;
    const double pop_a  = 4.5; // this is positive, the minus sign is already included in code
    
#warning dt and Rg filters could be run already while reading input files!
    
    const size_t min_selected = 1; // 2
    const size_t min_grains_per_pixel = 1; // 2
    
    ColDenDataContainer colden;
    
    for (int inputFileID=1; inputFileID<argc; ++inputFileID) {
        
        FILE * fp = fopen(argv[inputFileID],"r");
        if (fp) {
            ORSA_DEBUG("reading file [%s]",argv[inputFileID]);
        } else {
            ORSA_DEBUG("cannot open file [%s]",argv[inputFileID]);
            exit(0);
        }
        
        char line[4096];
        double r0,rF,AF,dt;
        double lon,lat;
        double pos_X,pos_Y,pos_Z;
        double pos_sun,pos_orbit_pole,pos_orbit_plane;
        double pos_orbit_velocity, /* pos_orbit_pole same as above */ pos_sunish;
        double pos_earth,pos_RA,pos_Dec;
#warning keep fields in sync with CometGrainsDynamics.cpp
        while (fgets(line,4096,fp)) {
            gmp_sscanf(line,
                       "%*s   %*s %*s   %lf %*s   %lf %lf   %lf %lf %lf   %*s   %lf %lf %lf   %lf %lf %lf    %lf %*s %lf    %lf %lf %lf   %*s",
                       &dt,
                       &lon,
                       &lat,
                       &r0,
                       &rF,
                       &AF,
                       //
                       &pos_X,
                       &pos_Y,
                       &pos_Z,
                       //
                       &pos_sun,
                       &pos_orbit_pole,
                       &pos_orbit_plane,
                       //
                       &pos_orbit_velocity,
                       // &pos_orbit_pole same as above
                       &pos_sunish,
                       //
                       &pos_earth,
                       &pos_RA,
                       &pos_Dec);
            dt = orsa::FromUnits(dt,orsa::Unit::DAY);
            r0 = orsa::FromUnits(r0,orsa::Unit::METER);
            rF = orsa::FromUnits(rF,orsa::Unit::METER);
            AF = orsa::FromUnits(AF,orsa::Unit::METER,2);
            //
            lon *= orsa::degToRad();
            lat *= orsa::degToRad();
            //
            pos_X = orsa::FromUnits(pos_X,orsa::Unit::KM);
            pos_Y = orsa::FromUnits(pos_Y,orsa::Unit::KM);
            pos_Z = orsa::FromUnits(pos_Z,orsa::Unit::KM);
            //
            pos_sun         = orsa::FromUnits(pos_sun,orsa::Unit::KM);
            pos_orbit_pole  = orsa::FromUnits(pos_orbit_pole,orsa::Unit::KM);
            pos_orbit_plane = orsa::FromUnits(pos_orbit_plane,orsa::Unit::KM);
            //
            pos_orbit_velocity = orsa::FromUnits(pos_orbit_velocity,orsa::Unit::KM);
            // pos_orbit_pole same as above
            pos_sunish = orsa::FromUnits(pos_sunish,orsa::Unit::KM);
            //
            pos_earth = orsa::FromUnits(pos_earth,orsa::Unit::KM);
            pos_RA    = orsa::FromUnits(pos_RA,orsa::Unit::KM);
            pos_Dec   = orsa::FromUnits(pos_Dec,orsa::Unit::KM);
            //
            if (0) {
                // test filter
                if (lat < -10.0*orsa::degToRad()) continue;
                if (lat > +10.0*orsa::degToRad()) continue;
                if (lon < -10.0*orsa::degToRad()) continue;
                if (lon > +10.0*orsa::degToRad()) continue;
            }
            //
            ColDenData data;
            //
            data.r0 = r0;
            data.rF = rF;
            data.AF = AF;
            data.dt = dt;
            //
            data.lon = lon;
            data.lat = lat;
            //
            data.pos_X = pos_X;
            data.pos_Y = pos_Y;
            data.pos_Z = pos_Z;
            //
            data.pos_sun = pos_sun;
            data.pos_orbit_pole = pos_orbit_pole;
            data.pos_orbit_plane = pos_orbit_plane;
            //
            data.pos_orbit_velocity = pos_orbit_velocity;
            // pos_orbit_pole same as above
            data.pos_sunish = pos_sunish;
            //
            data.pos_earth = pos_earth;
            data.pos_RA    = pos_RA;
            data.pos_Dec   = pos_Dec;
            //
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
    
    if (1) {
        
        // 2D
        typedef std::vector< std::vector< std::list<double> > > histo_2D_type;
        histo_2D_type histo_2D;
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
                        const double  pos_i = (*it).pos_X; // pos_orbit_velocity; // pos_RA; // pos_X; // pos_sun // UPDATE THIS
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
                        const double  pos_j = (*it).pos_Y; // .pos_sunish; // pos_Dec; // pos_Y; // pos_orbit_plane // UPDATE THIS
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
                    std::list<double>::const_iterator it = histo_2D[i][j].begin();
                    while (it != histo_2D[i][j].end()) {
                        val += (*it);
                        ++it;
                    }
                    if (val>0.0) {
                        ++populatedPixels;
                        fprintf(fp,"%3i %3i %.3e\n",i,j,val);
                    }
                }
            }
        }
        fclose(fp);
        
        ORSA_DEBUG("entries: %i   zero-area: %i   in-field: %i   pixels: %i (%.1f\%)",colden_size,zeroArea,inField,populatedPixels,100.0*(double)populatedPixels/(4*N*N));
        
    }
    
    if (0) {
        
        // 3D
        typedef std::vector< std::vector< std::vector< std::list<double> > > > histo_3D_type;
        histo_3D_type histo_3D;
        histo_3D.resize(2*N);
        for (size_t i=0; i<2*N; ++i) {
            histo_3D[i].resize(2*N);
            for (size_t j=0; j<2*N; ++j) {
                histo_3D[i][j].resize(2*N);
            }
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
                        const double  pos_i = (*it).pos_X; // pos_RA; // pos_X; // pos_sun // UPDATE THIS
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
                        const double  pos_j = (*it).pos_Y; // pos_Dec; // pos_Y; // pos_orbit_plane // UPDATE THIS
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
                        const double  pos_k = (*it).pos_Z; // pos_Dec; // pos_Y; // pos_orbit_plane // UPDATE THIS
                        const int k = N+pos_k/pixelScale;
                        if (k<0) {
                            ++it;
                            continue;
                        }
                        if (k>=2*N) {
                            ++it;
                            continue;
                        }
                        //
                        const double scaled_area = factor*(*it).AF;
                        
                        // choose one
                        //
                        histo_3D[i][j][k].push_back(scaled_area);
                        //
                        // rho (projected distance from nucleus) is for testing 1/rho
                        /* const double rho = sqrt(orsa::square(pos_i)+orsa::square(pos_j));
                           histo_3D[i][j][k].push_back(scaled_area*rho);
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
        FILE * fp = fopen("histo_3D.3D","w");
        fprintf(fp,"x y z area\n");
        for (int i=0; i<2*N; ++i) {
            for (int j=0; j<2*N; ++j) {
                for (int k=0; k<2*N; ++k) {
                    if (histo_3D[i][j][k].size()>=min_grains_per_pixel) {
                        double val = 0.0;
                        std::list<double>::const_iterator it = histo_3D[i][j][k].begin();
                        while (it != histo_3D[i][j][k].end()) {
                            val += (*it);
                            ++it;
                        }
                        if (val>0.0) {
                            ++populatedPixels;
                            // fprintf(fp,"%3i %3i %3i %.3e\n",i,j,k,val);
                            fprintf(fp,"%.f %.f %.f %.3e\n",
                                    orsa::FromUnits((i-N)*pixelScale,orsa::Unit::KM,-1),
                                    orsa::FromUnits((j-N)*pixelScale,orsa::Unit::KM,-1),
                                    orsa::FromUnits((k-N)*pixelScale,orsa::Unit::KM,-1),
                                    val);
                        }
                    }
                }
            }
        }
        fclose(fp);
        
        ORSA_DEBUG("entries: %i   zero-area: %i   in-field: %i   pixels: %i (%.1f\%)",colden_size,zeroArea,inField,populatedPixels,100.0*(double)populatedPixels/(8.0*N*N*N));
        
    }
    
    return 0;
}
