#ifndef _ORSA_SOLAR_SYSTEM_DATA_
#define _ORSA_SOLAR_SYSTEM_DATA_

namespace orsaSolarSystem {
  
    // as in orsa/unit.h, all units in MKS
  
    // data from: http://ssd.jpl.nasa.gov/?constants
    // updated from DE431 doc: https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
  
    class Data {
    public:
    
        /* enum ObjectID {
           MERCURY,
           VENUS,
           EARTH,
           MOON,
           EARTHMOON,
           MARS,
           JUPITER,
           SATURN,
           URANUS,
           NEPTUNE,
           PLUTO
           };
        */
    
        // GM coefficients
    
        inline static const double & GMSun() {
            static const double _GM = 1.32712440041939400e20;
            return _GM;
        }
    
        inline static const double & GMMercury() {
            static const double sunMassRatio = 6023682.155592; 
            static const double _GM = GMSun()/sunMassRatio;
            return _GM;
        }
    
        inline static const double & GMVenus() {
            static const double sunMassRatio = 408523.718658; 
            static const double _GM = GMSun()/sunMassRatio;
            return _GM;
        }
        
        /*
        inline static const double & EarthMoonMassRatio() {
            static const double EMratio = 81.30059;
            return EMratio;
        }
        */
        
        /*
        inline static const double & GMEarth() {
            static const double _GM = GMEarthMoon()*(EarthMoonMassRatio()/(EarthMoonMassRatio()+1));
            return _GM;
        }
        */
        
        inline static const double & GMEarth() {
            static const double sunMassRatio = 332946.048834; 
            static const double _GM = GMSun()/sunMassRatio;
            return _GM;
        }
        
        inline static const double & GMMoon() {
            static const double sunMassRatio = 27068703.241203; 
            static const double _GM = GMSun()/sunMassRatio;
            return _GM;
        }
    
        inline static const double & GMEarthMoon() {
            static const double sunMassRatio = 328900.559810; 
            static const double _GM = GMSun()/sunMassRatio;
            return _GM;
        }
    
        inline static const double & GMMars() {
            static const double sunMassRatio = 3098703.590291; 
            static const double _GM = GMSun()/sunMassRatio;
            return _GM;
        }
    
        inline static const double & GMJupiter() {
            static const double sunMassRatio = 1047.348625; 
            static const double _GM = GMSun()/sunMassRatio;
            return _GM;
        }
    
        inline static const double & GMSaturn() {
            static const double sunMassRatio = 3497.901768; 
            static const double _GM = GMSun()/sunMassRatio;
            return _GM;
        }
    
        inline static const double & GMUranus() {
            static const double sunMassRatio = 22902.981613; 
            static const double _GM = GMSun()/sunMassRatio;
            return _GM;
        }
    
        inline static const double & GMNeptune() {
            static const double sunMassRatio = 19412.259776; 
            static const double _GM = GMSun()/sunMassRatio;
            return _GM;
        }
    
        inline static const double & GMPluto() {
            static const double sunMassRatio = 1.35836683768617e8; 
            static const double _GM = GMSun()/sunMassRatio;
            return _GM;
        }
    
        // masses, derived from GM
    
        inline static const double & MSun() {
            static const double _M = GMSun()/orsa::Unit::G();
            return _M;
        }
    
        inline static const double & MMercury() {
            static const double _M = GMMercury()/orsa::Unit::G();
            return _M;
        }
    
        inline static const double & MVenus() {
            static const double _M = GMVenus()/orsa::Unit::G();
            return _M;
        }
    
        inline static const double & MEarth() {
            static const double _M = GMEarth()/orsa::Unit::G();
            return _M;
        }
    
        inline static const double & MMoon() {
            static const double _M = GMMoon()/orsa::Unit::G();
            return _M;
        }
    
        inline static const double & MEarthMoon() {
            static const double _M = GMEarthMoon()/orsa::Unit::G();
            return _M;
        }
    
        inline static const double & MMars() {
            static const double _M = GMMars()/orsa::Unit::G();
            return _M;
        }
    
        inline static const double & MJupiter() {
            static const double _M = GMJupiter()/orsa::Unit::G();
            return _M;
        }
    
        inline static const double & MSaturn() {
            static const double _M = GMSaturn()/orsa::Unit::G();
            return _M;
        }
    
        inline static const double & MUranus() {
            static const double _M = GMUranus()/orsa::Unit::G();
            return _M;
        }
    
        inline static const double & MNeptune() {
            static const double _M = GMNeptune()/orsa::Unit::G();
            return _M;
        }
    
        inline static const double & MPluto() {
            static const double _M = GMPluto()/orsa::Unit::G();
            return _M;
        }

        // other constants
    
        inline static const double & REarth() {
            static const double _R = 6378136.3;
            return _R;
        }
    
    };  
  
}; // namespace orsaSolarSystem

#endif 
