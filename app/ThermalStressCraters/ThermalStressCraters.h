#ifndef TSC_H
#define TSC_H

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

#include <orsa/debug.h>
#include <orsa/matrix.h>
#include <orsa/orbit.h>
#include <orsa/util.h>

#include <orsaSolarSystem/attitude.h>
#include <orsaSolarSystem/data.h>
#include <orsaSolarSystem/print.h>

#include <orsaSPICE/spice.h>
#include <orsaSPICE/spiceBodyTranslationalCallback.h>
#include <orsaSPICE/spiceBodyRotationalCallback.h>

// All thermal equations and notation from:
// Spencer, Lebofsky, Sykes 1989.
// "Systematic Biases in Radiometric Diameter Determinations"
// Icarus 78, 337-254.

// all in SI units

inline double bondAlbedo() { return 0.05; } // 0.05 // 0.20

inline double emissivity() { return 0.90; } // 0.90

inline double sigma() { return 5.67040e-8; } // Stefan-Boltzmann

inline double solar() { return 1370.0; } // W m^-2 at 1 AU

#warning update rotationPeriod
inline double rotationPeriod() { return 10.0*3600.0; } // s

#warning for omega, we need the sidereal or solar rotation period?
inline double omega() { return orsa::twopi()/rotationPeriod(); } // s^-1

// k
// const double thermalConductivity = 1.71e-4; // W m^-1 K^-1
inline double thermalConductivity() { return 1.90e-3; } // W m^-1 K^-1

// rho
inline double density() { return 1750.0; } // kg/m^3

// c
inline double specificHeatCapacity() { return 750.0; } // J kg^-1 K^-1

// Gamma
inline double thermalInertia() { return sqrt(thermalConductivity()*density()*specificHeatCapacity()); } // Joule m^-2 s^-1/2 K^-1

// l_s
inline double skinDepth() { return sqrt(thermalConductivity()/(density()*specificHeatCapacity()*omega())); } // m

// sub-solar Temperature, K
inline double Tss(const double & hdist) {
    return pow((1.0-bondAlbedo())*(solar()/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2))/(emissivity()*sigma()),0.25);
}

inline double theta(const double hdist) { return thermalInertia()*sqrt(omega())/(emissivity()*sigma()*pow(Tss(hdist),3)); } 

class Slice {
public:
    double T;
};

typedef std::vector<Slice> ThermalData;

typedef std::vector<ThermalData>  History;

class Profile {
public:
    Profile(const int numSlices_in) :
        numSlices(numSlices_in),
        ls(skinDepth()) {
        data.resize(numSlices+1);
    }
public:
    void step(const double dx,
              const double dt,
              const double Fs) {
    
        // const double dX = dx/skinDepth();
        const double dX = dx/ls;
        
        // ORSA_DEBUG("dX: %g",dX);
        
        const ThermalData old_data = data;
        
        for (int k=0; k<numSlices; ++k) {
            if (old_data[k].T < 0.0) {
                ORSA_DEBUG("problems: negative temperature...");
                exit(0);
            }
        }
        
        // ORSA_DEBUG("dX: %g dt: %g tI: %g Fs: %g Om: %g pw: %g old_data[0].T: %g",dX,dt,thermalInertia(),Fs,omega(),pow(old_data[0].T,4),old_data[0].T);
        
        // surface
        data[0].T =
            old_data[0].T
            + 2*dt*omega()/(dX*dX)*(old_data[1].T-old_data[0].T) 
            - 2*dt*sqrt(omega())/(thermalInertia()*dX)*(emissivity()*sigma()*pow(old_data[0].T,4)-(1.0-bondAlbedo())*Fs);
        // orsa::check(data[0].T);
        
        // explicitly skip k==0 and k==numSlices
        for (int k=1; k<numSlices; ++k) {
            data[k].T = 
                old_data[k].T
                + (dt*omega())/(dX*dX)*(old_data[k+1].T - 2*old_data[k].T + old_data[k-1].T);
            // orsa::check(data[k].T);
        }
        
    }
public:
    const int numSlices;
    const double ls; // skinDepth
public:
    ThermalData data;
protected:
    // ThermalData old_data;
};

bool ComputePeriodicThermalHistory(History & history,
                                   const size_t numSlices, // typically equal to numSkinDepths*numSlicesPerSkinDepth
                                   const double & initial_deep_T, // K
                                   const std::vector<double> & Fs,
                                   const double dx,
                                   const double dt,
                                   const unsigned int & history_skip=1,
                                   const double & stability_eps=1.0e-4,
                                   const double & convergence_eps=1.0e-9) {
    
    // Profile profile(numSkinDepths*numSlicesPerSkinDepth);  
    Profile profile(numSlices);  
  
    // History history;
    // const unsigned int history_skip = numTimeStepsPerRotation/10000; // save memory, same results, affects print-out density of results too
    // const unsigned int history_skip = 1; // save memory, similar results, affects print-out density of results too
    // ORSA_DEBUG("history_skip: %i",history_skip);
    history.clear();
    
    // initial temperature at max depth
    // double deep_T = 150.0; // K
    // double deep_T = 0.0; // K
    double deep_T = initial_deep_T;
    
    for (unsigned int k=0; k<profile.data.size(); ++k) {
        profile.data[k].T = deep_T;
    }
    
    bool converged=false;
    while (!converged) {
        
        ORSA_DEBUG("deep_T: %g [K]",deep_T);
        
        ORSA_DEBUG("dt: %f [s]   (with history_skip=%i then saved dt: %f [s])   dx: %g [m]",dt,history_skip,history_skip*dt,dx);
        
        bool stable=false;
        History old_history;
        unsigned int cycles=0;
        // const unsigned int stable_max_try=32;
        while (!stable) {
            old_history = history;
            history.clear();
            /* for (int p=0; p<=numTimeStepsPerRotation*totalRotations; ++p) {
               const orsa::Time t = t0+p*dt;
               // orsaSolarSystem::print(t);
               orsa::Vector rVesta;
               bg->getInterpolatedPosition(rVesta,vesta.get(),t);
               orsa::Vector rSun;
               bg->getInterpolatedPosition(rSun,sun.get(),t);
               const orsa::Vector dr = (rVesta-rSun);
               const orsa::Vector u_surface_to_sun = (-dr).normalized();
               const double hdist = dr.length();
               const double scaled_solar = solar/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2);
               const orsa::Matrix l2g = orsa::localToGlobal(vesta.get(),bg.get(),t);
               const orsa::Vector u_global_normal_surface_element = l2g*normal_local;
               const double Fs = scaled_solar*std::max(0.0,u_surface_to_sun*u_global_normal_surface_element);
               // ORSA_DEBUG("Fs: %g   scaled_solar: %g   p: %i   dt: %g",Fs,scaled_solar,p,dt.get_d());
               profile.step(dx,dt.get_d(),Fs);
               if (p%history_skip==0) history.push_back(profile.data);	
               }
            */
            for (size_t p=0; p<Fs.size(); ++p) {
                /*
                   if (p==0) {
                   // search stable first profile
                   size_t iter=0;
                   while (1) {
                   ++iter;
                   Profile old_profile = profile;
                   profile.step(dx,dt,Fs[p]);
                   stable=true;
                   for (unsigned int k=0; k<profile.data.size(); ++k) {
                   if (!stable) break;
                   if (fabs((profile.data[k].T-old_profile.data[k].T)/old_profile.data[k].T) > stability_eps) {
                   stable=false;
                   break;
                   }
                   }
                   if (stable) {
                   ORSA_DEBUG("first step stable after %i cycles",iter);
                   break;
                   }
                   }
                   }
                */
                
                profile.step(dx,dt,Fs[p]);
                if (p%history_skip==0) {
                    history.push_back(profile.data);
                }
            }
            
            // const double stable_eps = 1.0e-4;
            // const double stable_eps = 1.0e-2;
            
            stable=true;
            if (old_history.size()==0) {
                ORSA_DEBUG("stable search did not converge yet, cycle # %03i",
                           cycles);
                stable=false;
            } else {
                for (unsigned int k=0; k<history[0].size(); ++k) {
                    if (!stable) break;
                    for (unsigned int j=0; j<history.size(); ++j) {
                        if (!stable) break;
                        if (fabs((history[j][k].T-old_history[j][k].T)/old_history[j][k].T) > stability_eps) {
                            ORSA_DEBUG("stable search did not converge yet, cycle # %03i   delta = %.3e > %.3e at bin: %i/%i time: %i/%i",
                                       cycles,
                                       fabs((history[j][k].T-old_history[j][k].T)/(orsa::epsilon()+old_history[j][k].T)),
                                       stability_eps,
                                       k,
                                       history[0].size(),
                                       history_skip*j,
                                       history_skip*history.size());
                            stable=false;
                            break;
                        }   
                    }
                }
            }
            ++cycles;
            // if (cycles >= stable_max_try) break;
      
            // when debugging
            // #warning REMOVE!
            // if (cycles==10) exit(0);
        }
        //
        if (stable) ORSA_DEBUG("stable after %i cycles",cycles);
    
        // determine constant value of deep temperature
        std::vector<double> average;
        average.resize(history[0].size());
        for (unsigned int k=0; k<history[0].size(); ++k) {
            average[k]=0;
            for (unsigned int j=0; j<history.size(); ++j) {
                average[k] += history[j][k].T;
            }
            average[k] /= history.size();
            // ORSA_DEBUG("average[%03i] = %.6f",k,average[k]);
        }
    
        if (stable) {
            const double converged_eps = 1.0e-2;
            // const double converged_eps = 1.0e-1;
            converged=true;
            for (unsigned int k=0; k<average.size(); ++k) {
                if (fabs((average[k]-average[0])/average[0]) > converged_eps) {
                    ORSA_DEBUG("not converged: %g > %g, k=%i",fabs((average[k]-average[0])/average[0]),converged_eps,k);
                    converged=false;
                    break;
                }
            }
        } else {
            converged=false;
        }
        
        // base filename
        // char filename[4096];
        
        // physical output
        /* if (1) {
           static unsigned int fileCounter=0;
           ++fileCounter;
           sprintf(filename,"temperature_orbital_%+05.2f_%g_%sgcm2_v%03i.dat",lat*orsa::radToDeg(),lon*orsa::radToDeg(),argv[3],fileCounter);
           ORSA_DEBUG("writing output file %s",filename);
           FILE * fp = fopen(filename,"w");
           for (unsigned int j=0; j<history.size(); ++j) {
           const orsa::Time t = t0+j*history_skip*dt;
           orsa::Vector rVesta;
           bg->getInterpolatedPosition(rVesta,vesta.get(),t);
           orsa::Vector rSun;
           bg->getInterpolatedPosition(rSun,sun.get(),t);
           const orsa::Vector dr = (rVesta-rSun);
           const orsa::Vector u_surface_to_sun = (-dr).normalized();
           const double hdist = dr.length();
           const double scaled_solar = solar/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2);
           const orsa::Matrix l2g = orsa::localToGlobal(vesta.get(),bg.get(),t);
           const orsa::Matrix g2l = orsa::globalToLocal(vesta.get(),bg.get(),t);
           const orsa::Vector u_global_normal_surface_element = l2g*normal_local;
           // b2s vector to compute sub-solar longitude
           orsa::Vector b2s_local = g2l*(rSun-rVesta);
           const double subSolarLongitude = atan2(b2s_local.getY(),b2s_local.getX());
           const double Fs = scaled_solar*std::max(0.0,u_surface_to_sun*u_global_normal_surface_element);
           // GRaND coefficient: weighted average of temperature with depth, weighting function negative exponential
           double GRaND_coefficient = 0.0;
           {
           double sum_weight = 0.0;
           for (unsigned int k=0; k<history[j].size(); ++k) {
           GRaND_coefficient += exp(-(k*dx/d0))*history[j][k].T;
           sum_weight += exp(-(k*dx/d0));
           }
           GRaND_coefficient /= sum_weight;
           // ORSA_DEBUG("GRaND: %g",GRaND_coefficient);
           }
           //
           fprintf(fp,"%10.5f %7.3f %10.6f %10.6f %10.6f %10.3f %10.6f %.2f %.2f\n",
           orsaSolarSystem::timeToJulian(t), // j*history_skip*dt.get_d()/rotationPeriod,
           orsa::radToDeg()*fmod((lon+orsa::pi()-subSolarLongitude)+2*orsa::twopi(),orsa::twopi()), // centered at noon
           history[j][0].T,
           history[j][history[j].size()-2].T, // -2 because -1 is never changed by thermal algo
           orsa::FromUnits(hdist,orsa::Unit::AU,-1),
           Fs,
           GRaND_coefficient,
           lat*orsa::radToDeg(),
           lon*orsa::radToDeg());
           }
           fclose(fp);
           //
           char cmd[1024];
           sprintf(cmd,"cp %s temperature_orbital_%+05.2f_%g_%sgcm2_latest.dat",filename,lat*orsa::radToDeg(),lon*orsa::radToDeg(),argv[3]);
           ORSA_DEBUG("executing: [%s]",cmd);
           int retval = system(cmd);
           if (retval != 0) ORSA_DEBUG("problems with the system call...");
           }
        */
        
        // correction
        if (!converged) {
            const double new_deep_T = average[0];
            if (fabs((new_deep_T-deep_T)/deep_T) < convergence_eps) {
                ORSA_DEBUG("temperature correction not converging, try to change integration steps, exiting");
                exit(0);
            }
            ORSA_DEBUG("new initial temperature: %.3f",new_deep_T);
            for (unsigned int l=0; l<profile.data.size(); ++l) {
                profile.data[l].T = new_deep_T;
            }
            deep_T = new_deep_T;
        } else {
            ORSA_DEBUG("converged");
            /* char cmd[1024];
               sprintf(cmd,"cp %s temperature_orbital_%+05.2f_%g_%sgcm2_final.dat",filename,lat*orsa::radToDeg(),lon*orsa::radToDeg(),argv[3]);
               ORSA_DEBUG("executing: [%s]",cmd);
               int retval = system(cmd);
               if (retval != 0) ORSA_DEBUG("problems with the system call...");
               ORSA_DEBUG("converged");
            */
        }
        
    }
    
    return converged;
}

bool CraterShape(double & h, /* elevation, from 0.0 (rim) to -d (center) */
                 double & dhdr, /* slope at the point r */
                 const double & r, /* point in crater, from 0 (center) to R=D/2 (rim=edge) */
                 const double & D, /* Diameter */
                 const double & d, /* depth */
                 const double & alpha0, /* slope at crater center = tan(slope_angle)*/
                 const double & alphaR /* slope at crater rim */) {
    
    const double R = 0.5*D;
    
    if (alpha0<0.0) {
        ORSA_DEBUG("negative slope at center of the crater");
        return false;
    }
    
    if (alpha0>=(d/R)) {
        ORSA_DEBUG("too steep at center of the crater");
        return false;
    }
    
    if (alphaR<=(d/R)) {
        ORSA_DEBUG("too shallow at the rim of the crater");
        return false;
    }
    
    if (alphaR>=orsa::pi()) {
        ORSA_DEBUG("too steep at the rim of the crater");
        return false;
    }

    if (r<0.0) {
        ORSA_DEBUG("problem: negative radius");
        return false;
    }

    if (r>=R) {
        ORSA_DEBUG("point is outside crater");
        h = 0.0;
        dhdr = 0.0;
        return true;
    }
    
    const double gamma = (alphaR-d/R)/(d/R-alpha0);
    h = alpha0*r+(alphaR-alpha0)*pow(r/R,gamma+1)*R/(gamma+1)-d;
    dhdr = alpha0+(alphaR-alpha0)*pow(r/R,gamma);
    
    return true;
}

double SolarDiskFraction(const double & heliocentricDistance,
                         const double & solarCenterElevation /* (signed) angle from sun disk center to horizon */) {

    // test: zero solar size
    /* if (solarCenterElevation>0.0) {
       return 1.0;
       } else {
       return 0.0;
       }
    */
    
    const double solarRadius = orsa::FromUnits(6.955e8,orsa::Unit::METER);
    const double alpha = asin(solarRadius/heliocentricDistance); // solar radius angle
    if (solarCenterElevation>alpha) {
        return 1.0;
    } else if (solarCenterElevation<-alpha) {
        return 0.0;
    } else {
        const double var = solarCenterElevation/alpha;
        return (0.5 + (var*sqrt(1.0-var*var) + asin(var))/orsa::pi());
    }
}


#endif // TSC_H
