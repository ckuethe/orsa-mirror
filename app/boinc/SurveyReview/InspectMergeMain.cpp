#include "grain.h"

#include <orsa/debug.h>
#include <orsa/print.h>
#include <orsa/statistic.h>

#include <orsaSolarSystem/datetime.h>
#include <orsaSolarSystem/print.h>

#include "binomial.h"
#include "grain.h"
#include "SurveyReview.h"
#include "skycoverage.h"
#include "eta.h"
#include "fit.h"

// SQLite3
#include "sqlite3.h"

osg::ref_ptr<SkyCoverage> skyCoverage;

osg::ref_ptr<orsa::RNG> rnd;

osg::ref_ptr<orsa::BodyGroup> bg;
osg::ref_ptr<orsa::Body> sun;
osg::ref_ptr<orsa::Body> earth;
osg::ref_ptr<orsa::Body> moon; 
orsa::Orbit earthOrbit;

osg::ref_ptr<LinearVar> var_ra;
osg::ref_ptr<LinearVar> var_dec;

osg::ref_ptr<LinearVar> var_xy; // same var used twice for xy (ecliptic plane)

// epoch for sky projection
orsa::Cache<orsa::Time> epoch;

class PlotStatsElement_WS : public orsa::WeightedStatistic<double> {
public:
    void insert(const double & x) {
        orsa::WeightedStatistic<double>::insert(x,1.0);
    }
};

class PlotStatsElement_Sum : public osg::Referenced {
public:
    PlotStatsElement_Sum() : osg::Referenced(), _s(0.0), _n(0) { }
protected:
    virtual ~PlotStatsElement_Sum() { }
public:
    void insert(const double & x) {
        _s += x;
        ++_n;
    }
public:
    const double & sum() const {
        return _s;
    }
public:
    const mpz_class & entries() const {
        return _n;
    }
protected:
    double _s;
    mpz_class _n;
};

template <typename E, typename B> class PlotStats : public BinStats<E> {
public:
    PlotStats(const std::vector< osg::ref_ptr<Var> > & varDefinition) :
        BinStats<E>(varDefinition) { }
public:
    bool insert(const std::vector<double> & xVector,
                const double & val) {
        if (xVector.size() != B::var.size()) {
            ORSA_DEBUG("dimension mismatch");
            return false;
        }
        std::vector<size_t> binVector;
        if (!B::bin(binVector,xVector)) {
            return false;
        }
        const mpz_class idx = B::index(binVector);
        if (B::data[idx].get()==0) {
            // lazy allocation
            B::data[idx] = new E;
        }
        B::data[idx]->insert(val);
        // ORSA_DEBUG("xV: %g %g bV: %i %i",xVector[0],xVector[1],binVector[0],binVector[1]);
        return true;         
    }
};

typedef PlotStats<PlotStatsElement_WS,  BinStats<PlotStatsElement_WS>  > PlotStats_WS;
typedef PlotStats<PlotStatsElement_Sum, BinStats<PlotStatsElement_Sum> > PlotStats_Sum;

#warning UPDATE writeOutputFile... code with faster version in app/dawnorbit/postproc.cpp

// Nsub is the number of sub-bins in each x,y bin -- to account for missing zero entries
// i.e. if writing a,e then Nsub=Ni*Nnode*Nperi*Nm=18*12*12*12...
void writeOutputFile(const std::string & filename,
                     PlotStats_WS * plotStats,
                     const LinearVar * var_x,
                     const LinearVar * var_y,
                     const unsigned int Nsub) {
    
    FILE * fp = fopen(filename.c_str(),"w");    
    
    std::vector<double> xVector;
    xVector.resize(2);
    
    for (unsigned j=0; j<var_x->size(); ++j) {
        for (unsigned k=0; k<var_y->size(); ++k) {
            
            xVector[0] = var_x->start+var_x->incr*(j+0.5);
            xVector[1] = var_y->start+var_y->incr*(k+0.5);
            
            std::vector<size_t> binVector;
            if (plotStats->bin(binVector,xVector)) {
                const PlotStatsElement_WS * e = plotStats->stats(plotStats->index(binVector));
                if (e) {
                    gmp_fprintf(fp,
                                "%8g %8g %8.6f %8.6f %6Zi\n",
                                xVector[0],
                                xVector[1],
                                e->average(),
                                e->average()*e->entries().get_d()/Nsub, // divide by Nsub to account for missing zero entries
                                e->entries().get_mpz_t());
#warning the two should be the same for a complete analysis, because over several years, each a,e bin has been observed at least once!
                }
            }
        }
    }
    
    fclose(fp);
}

// Nsub is the number of sub-bins in each x,y bin -- to account for missing zero entries
// i.e. if writing a,e then Nsub=Ni*Nnode*Nperi*Nm=18*12*12*12...
void writeOutputFile(const std::string & filename,
                     PlotStats_Sum * plotStats,
                     const LinearVar * var_x,
                     const LinearVar * var_y,
                     const unsigned int Nsub,
                     const bool sky_norm=false) {
    
    FILE * fp = fopen(filename.c_str(),"w");    
    
    std::vector<double> xVector;
    xVector.resize(2);
    
    double area=1.0;
    
    for (unsigned j=0; j<var_x->size(); ++j) {
        for (unsigned k=0; k<var_y->size(); ++k) {
            
            xVector[0] = var_x->start+var_x->incr*(j+0.5);
            xVector[1] = var_y->start+var_y->incr*(k+0.5);

#warning bin-by-bin normalization (divide by area), find better way to implement
            if (sky_norm) {
                area = fabs(var_x->incr*15.0*orsa::degToRad()*
                            (sin(orsa::degToRad()*(var_y->start+var_y->incr*(k))) -
                             sin(orsa::degToRad()*(var_y->start+var_y->incr*(k+1)))));
            }
            
            // ORSA_DEBUG("xV: %g %g  area: %g",xVector[0],xVector[1],area);
            
            std::vector<size_t> binVector;
            if (plotStats->bin(binVector,xVector)) {
                const PlotStatsElement_Sum * e = plotStats->stats(plotStats->index(binVector));
                if (e) {
                    gmp_fprintf(fp,
                                "%8g %8g %8.6f %8.6f %6Zi\n",
                                xVector[0],
                                xVector[1],
                                (e->sum()/area),
                                (e->sum()/area)*e->entries().get_d()/Nsub, // divide by Nsub to account for missing zero entries
                                e->entries().get_mpz_t());
#warning the two should be the same for a complete analysis, because over several years, each a,e bin has been observed at least once!
                }
            }
        }
    }
    
    fclose(fp);
}

// global vars, for use in the callbacks

// H16
osg::ref_ptr< PlotStats_WS > plotStats_ae_NEO_H16;
osg::ref_ptr< PlotStats_WS > plotStats_ae_PHO_H16;
//
osg::ref_ptr< PlotStats_WS > plotStats_ai_NEO_H16;
osg::ref_ptr< PlotStats_WS > plotStats_ai_PHO_H16;
//
osg::ref_ptr< PlotStats_WS > plotStats_aL_NEO_H16;
osg::ref_ptr< PlotStats_WS > plotStats_aL_PHO_H16;
//
osg::ref_ptr< PlotStats_Sum > plotStats_sky_NEO_H16;
osg::ref_ptr< PlotStats_Sum > plotStats_sky_PHO_H16;
//
osg::ref_ptr< PlotStats_Sum > plotStats_xy_NEO_H16;
osg::ref_ptr< PlotStats_Sum > plotStats_xy_PHO_H16;

// H18
osg::ref_ptr< PlotStats_WS > plotStats_ae_NEO_H18;
osg::ref_ptr< PlotStats_WS > plotStats_ae_PHO_H18;
//
osg::ref_ptr< PlotStats_WS > plotStats_ai_NEO_H18;
osg::ref_ptr< PlotStats_WS > plotStats_ai_PHO_H18;
//
osg::ref_ptr< PlotStats_WS > plotStats_aL_NEO_H18;
osg::ref_ptr< PlotStats_WS > plotStats_aL_PHO_H18;
//
osg::ref_ptr< PlotStats_Sum > plotStats_sky_NEO_H18;
osg::ref_ptr< PlotStats_Sum > plotStats_sky_PHO_H18;
//
osg::ref_ptr< PlotStats_Sum > plotStats_xy_NEO_H18;
osg::ref_ptr< PlotStats_Sum > plotStats_xy_PHO_H18;

int inspectCallback(void  * /* unused */,
                    int     /* ncols  */,
                    char ** col,
                    char ** /* colName */) {
    
    static mpz_class callID=0;
    
    const int z_a_min          = atoi(col[0]);
    const int z_a_max          = atoi(col[1]);
    const int z_e_min          = atoi(col[2]);
    const int z_e_max          = atoi(col[3]);
    const int z_i_min          = atoi(col[4]);
    const int z_i_max          = atoi(col[5]);
    const int z_node_min       = atoi(col[6]);
    const int z_node_max       = atoi(col[7]);
    const int z_peri_min       = atoi(col[8]);
    const int z_peri_max       = atoi(col[9]);
    const int z_M_min          = atoi(col[10]);
    const int z_M_max          = atoi(col[11]);
    //
    const int z_H              = atoi(col[12]);
    //
    const int N_NEO            = atoi(col[13]);
    const int N_PHO            = atoi(col[14]);
    const int NEO_in_field     = atoi(col[15]);
    const int PHO_in_field     = atoi(col[16]);
    //
    const double eta_NEO       = atof(col[17]);
    // const double sigma_eta_NEO = atof(col[18]);
    const double eta_PHO       = atof(col[19]);
    // const double sigma_eta_PHO = atof(col[20]);
    
    if ( (NEO_in_field < 0) ||
         (PHO_in_field < 0) ) {
        ORSA_DEBUG("skipping negative entry... is this a merged db?");
        return 0;
    }
    
#warning keep this filter updated
    if ( (z_H != 160) && (z_H != 180) ) {
        // quick exit, since we're saving only these two for now...
        return 0;
    }
    
    // local center bin
    const double center_a = 0.5*(z_a_max+z_a_min)*grain_a_AU;
    const double center_e = 0.5*(z_e_max+z_e_min)*grain_e;
    const double center_i = 0.5*(z_i_max+z_i_min)*grain_i_DEG;

    // size=2 should work for most cases
    std::vector<double> xVector; xVector.resize(2);
    // xVector.resize(varDefinition.size());
    
    if (z_H == 160) {
        xVector[0] = center_a;
        xVector[1] = center_e;
        plotStats_ae_NEO_H16->insert(xVector, eta_NEO);
        plotStats_ae_PHO_H16->insert(xVector, eta_PHO);
        //
        xVector[0] = center_a;
        xVector[1] = center_i;
        plotStats_ai_NEO_H16->insert(xVector, eta_NEO);
        plotStats_ai_PHO_H16->insert(xVector, eta_PHO);
    }
    
    if (z_H == 180) {
        xVector[0] = center_a;
        xVector[1] = center_e;
        plotStats_ae_NEO_H18->insert(xVector, eta_NEO);
        plotStats_ae_PHO_H18->insert(xVector, eta_PHO);
        //
        xVector[0] = center_a;
        xVector[1] = center_i;
        plotStats_ai_NEO_H18->insert(xVector, eta_NEO);
        plotStats_ai_PHO_H18->insert(xVector, eta_PHO);
    }
    
    // L
    std::vector<int> z_L_vec;
    {
        // assume step is the same for node,peri,M
        // #warning this is grain-size dependent!
        z_L_vec.push_back(((z_node_min+z_peri_min+z_M_min)%360000)/30000);
        z_L_vec.push_back(((z_node_max+z_peri_min+z_M_min)%360000)/30000);
        z_L_vec.push_back(((z_node_max+z_peri_max+z_M_min)%360000)/30000);
    }
    //
    for (unsigned int l=0; l<z_L_vec.size(); ++l) {
        const int index_L = z_L_vec[l]; // 0-11
        
        // assuming grain size for node=peri=M=L
        const double center_L = (0.5+index_L)*30.0;
        
        if (z_H == 160) {
            xVector[0] = center_a;
            xVector[1] = center_L;
            plotStats_aL_NEO_H16->insert(xVector, eta_NEO);
            plotStats_aL_PHO_H16->insert(xVector, eta_PHO);
        }
        
        if (z_H == 180) {
            xVector[0] = center_a;
            xVector[1] = center_L;
            plotStats_aL_NEO_H18->insert(xVector, eta_NEO);
            plotStats_aL_PHO_H18->insert(xVector, eta_PHO);
        }
    }
    
    // sky + xy
    {
        // don't use all binomial code for now, just the average value
        //
        // no more than 1.0 for now...
        // notice the "-NEO_in_field" at the end: missing pop = estimated pop - known pop
        /* double zpop=1.0;
           if ((eta_NEO > 0.01) && (eta_NEO <= 1.0)) {
           zpop=(((NEO_in_field+1)/eta_NEO)-1)-NEO_in_field; 
           }
           const unsigned int N_NEO_missing_pop_mean = lrint(std::min(1.0,zpop));    
           #warning TODO: sample 100 objects from pop bin, and project them on the sky, each one with weight N_NEO_missing_pop_mean/100
           #warning also weight by eta_detection, computed for each object one by one
        */
        const double prob_one_more = probabilityToFindOneMore(eta_NEO,NEO_in_field,true);
        
        const orsa::Time apparentMotion_dt_T = orsa::Time(0,0,1,0,0);
        //
        const double apparentMotion_dt = apparentMotion_dt_T.get_d();
        
        orsa::Vector observerPosition_epoch;
        orsa::Vector observerPosition_epoch_plus_dt;
        orsa::Vector sunPosition_epoch;
        orsa::Vector sunPosition_epoch_plus_dt;
        
        // using geocentric observer, so obs-position = earth positon 
        /* obsPosCB->getPosition(observerPosition_epoch,
           skyCoverage->obscode,
           skyCoverage->epoch);
           obsPosCB->getPosition(observerPosition_epoch_plus_dt,
           skyCoverage->obscode,
           skyCoverage->epoch+apparentMotion_dt_T);
        */

        /* bg->getInterpolatedPosition(observerPosition_epoch,
           earth.get(),
           skyCoverage->epoch);
           
           bg->getInterpolatedPosition(observerPosition_epoch_plus_dt,
           earth.get(),
           skyCoverage->epoch+apparentMotion_dt_T);
           
           bg->getInterpolatedPosition(sunPosition_epoch,
           sun.get(),
           skyCoverage->epoch);
           
           bg->getInterpolatedPosition(sunPosition_epoch_plus_dt,
           sun.get(),
           skyCoverage->epoch+apparentMotion_dt_T);
        */
        
        // computed at skyCoverage->epoch
        /* const orsa::Vector observerPosition_sk_epoch         = observerPosition_epoch;
           const orsa::Vector observerPosition_sk_epoch_plus_dt = observerPosition_epoch_plus_dt;
           const orsa::Vector sunPosition_sk_epoch              = sunPosition_epoch;
           const orsa::Vector sunPosition_sk_epoch_plus_dt      = sunPosition_epoch_plus_dt;
        */
        
        osg::ref_ptr<OrbitFactory> orbitFactory =
            new OrbitFactory(grain_a_AU*z_a_min,
                             grain_a_AU*z_a_max,
                             grain_e*z_e_min,
                             grain_e*z_e_max,
                             grain_i_DEG*z_i_min,
                             grain_i_DEG*z_i_max,
                             grain_node_DEG*z_node_min,
                             grain_node_DEG*z_node_max,
                             grain_peri_DEG*z_peri_min,
                             grain_peri_DEG*z_peri_max,
                             grain_M_DEG*z_M_min,
                             grain_M_DEG*z_M_max,
                             // grain_H* z_H,
                             // grain_H*(z_H+z_H_delta),
                             rnd.get(),
                             earthOrbit);
        
        unsigned int sky_N_NEO=0;
        unsigned int sky_N_PHO=0;
        
#warning this JD depends on the run, shoud be fixed in one file
        // JD from SRMJS...
#warning check for timescale info?
        const double JD = 2455650; // epoch of orbits
        const orsa::Time orbitEpoch = orsaSolarSystem::julianToTime(JD); 
        // ORSA_DEBUG("orbit epoch [below]");
        // orsaSolarSystem::print(orbitEpoch);
        
        bg->getInterpolatedPosition(observerPosition_epoch,
                                    earth.get(),
                                    epoch);
        
        bg->getInterpolatedPosition(observerPosition_epoch_plus_dt,
                                    earth.get(),
                                    epoch+apparentMotion_dt_T);
        
        bg->getInterpolatedPosition(sunPosition_epoch,
                                    sun.get(),
                                    epoch);
        
        bg->getInterpolatedPosition(sunPosition_epoch_plus_dt,
                                    sun.get(),
                                    epoch+apparentMotion_dt_T);
        
        orsa::Vector earthPosition_epoch;
        bg->getInterpolatedPosition(earthPosition_epoch,
                                    earth.get(),
                                    epoch);
        
        /* ORSA_DEBUG("earth position: %g %g %g",
           orsa::FromUnits(earthPosition_epoch.getX(),orsa::Unit::AU,-1),
           orsa::FromUnits(earthPosition_epoch.getY(),orsa::Unit::AU,-1),
           orsa::FromUnits(earthPosition_epoch.getZ(),orsa::Unit::AU,-1));
        */
        
        orsa::Vector moonPosition_epoch;
        bg->getInterpolatedPosition(moonPosition_epoch,
                                    moon.get(),
                                    epoch);
        
        // earth north pole
        const orsa::Vector northPole = (orsaSolarSystem::equatorialToEcliptic()*orsa::Vector(0,0,1)).normalized();

        const unsigned int target_num=16;
        for (unsigned int j=0; j<100; ++j) {
            // while (1) {
            
            osg::ref_ptr<OrbitID> orbit = orbitFactory->sample();
			
            if (!orbit->isNEO()) {
                continue;
            }
            ++sky_N_NEO;
            
            const bool isPHO = orbit->isPHO();
            if (isPHO) {
                ++sky_N_PHO;
            }
            
            const double orbitPeriod = orbit->period();
            
            /* obsPosCB->getPosition(observerPosition_epoch,
               skyCoverage->obscode,
               epoch);
               
               obsPosCB->getPosition(observerPosition_epoch_plus_dt,
               skyCoverage->obscode,
               epoch+apparentMotion_dt_T);
            */

            orsa::Vector r;
            
            const double original_M  = orbit->M;
            // 
            orbit->M = original_M + fmod(orsa::twopi() * (epoch-orbitEpoch).get_d() / orbitPeriod, orsa::twopi());
            orbit->relativePosition(r);
            orsa::Vector orbitPosition_epoch = r + sunPosition_epoch;
            //
            orbit->M = original_M + fmod(orsa::twopi() * (epoch+apparentMotion_dt_T-orbitEpoch).get_d() / orbitPeriod, orsa::twopi());
            orbit->relativePosition(r);
            orsa::Vector orbitPosition_epoch_plus_dt = r + sunPosition_epoch_plus_dt;
            //
            orbit->M = original_M;
            
            orsa::Vector dr_epoch          = (orbitPosition_epoch         - observerPosition_epoch);
            orsa::Vector dr_epoch_plus_dt  = (orbitPosition_epoch_plus_dt - observerPosition_epoch_plus_dt);
            
            // ++inField;
            
            const orsa::Vector orb2obs    = observerPosition_epoch - orbitPosition_epoch;
            const orsa::Vector obs2orb    = -orb2obs;
            const orsa::Vector orb2sun    = sunPosition_epoch      - orbitPosition_epoch;
            const orsa::Vector obs2sun    = sunPosition_epoch      - observerPosition_epoch;
            const orsa::Vector obs2moon   = moonPosition_epoch     - observerPosition_epoch;
            const double       phaseAngle = acos((orb2obs.normalized())*(orb2sun.normalized()));
            
            const orsa::Vector moon2obs = observerPosition_epoch - moonPosition_epoch;
            const orsa::Vector moon2sun =      sunPosition_epoch - moonPosition_epoch;
            
            // apparent magnitude V moved later, in H loop
            
            // apparent velocity
            const double U = acos(dr_epoch_plus_dt.normalized()*dr_epoch.normalized())/apparentMotion_dt;
            
            // airmass
            //
            // aliases
            // const orsa::Vector & earthPosition = earthPosition_epoch;
            // const orsa::Vector &   obsPosition = observerPosition_epoch;
            // const orsa::Vector     zenith = (obsPosition - earthPosition).normalized();
            // const orsa::Vector localEast  = orsa::externalProduct(northPole,zenith).normalized();
            //  const orsa::Vector localNorth = orsa::externalProduct(zenith,localEast).normalized();
            // const double obs2orb_zenith     =     zenith*obs2orb.normalized();
            // const double obs2orb_localEast  =  localEast*obs2orb.normalized();
            // const double obs2orb_localNorth = localNorth*obs2orb.normalized();
            // const double zenithAngle = acos(obs2orb_zenith);
            // const double airMass = ((observed||epochFromField)&&(zenithAngle<orsa::halfpi())?(1.0/cos(zenithAngle)):-1.0);
            // const double azimuth = fmod(orsa::twopi()+atan2(obs2orb_localEast,obs2orb_localNorth),orsa::twopi());
            // const double AM = ((epochFromField)&&(zenithAngle<orsa::halfpi())?(1.0/cos(zenithAngle)):100.0);
#warning NOTE: using airmass=1.0 because using generic geocentric observatory for now
            const double AM = 1.0;
            
            // const double solarAltitude = orsa::halfpi()-acos(zenith*obs2sun.normalized());
            // const double lunarAltitude = orsa::halfpi()-acos(zenith*obs2moon.normalized());
            // const double lunarElongation = acos(obs2moon.normalized()*obs2orb.normalized());
            // const double lunarPhase = acos(moon2obs.normalized()*moon2sun.normalized());
            //
            // const double SA = solarAltitude;
            // const double LA = lunarAltitude;
            // const double LE = lunarElongation;
            // const double LP = lunarPhase;
            
            // galactic latitude
            const orsa::Vector obs2orb_Equatorial = orsaSolarSystem::eclipticToEquatorial()*obs2orb;
            const orsa::Vector dr_equatorial = obs2orb_Equatorial.normalized();
            const double  ra = fmod(atan2(dr_equatorial.getY(),dr_equatorial.getX())+orsa::twopi(),orsa::twopi());
            const double dec = asin(dr_equatorial.getZ()/dr_equatorial.length());
            double l,b;
            orsaSolarSystem::equatorialToGalactic(l,b,ra,dec);
            // format longitude between -180 and +180 deg
            l = fmod(l+2*orsa::twopi(),orsa::twopi());
            if (l > orsa::pi()) l -= orsa::twopi();
            // const double galacticLongitude = l;
            // const double galacticLatitude  = b;
            const double GB = b;
            const double GL = l;
            // ecliptic coordinates
            /* const orsa::Vector dr = obs2orb.normalized();
               const double phi      = fmod(atan2(dr.getY(),dr.getX())+orsa::twopi(),orsa::twopi());
               const double theta    = asin(dr.getZ()/dr.length());
               const orsa::Vector dr_sun = obs2sun.normalized();
               const double phi_sun      = fmod(atan2(dr_sun.getY(),dr_sun.getX())+orsa::twopi(),orsa::twopi());
               const double theta_sun    = asin(dr_sun.getZ()/dr_sun.length());
               const double tmp_eclipticLongitude = fmod(phi-phi_sun+orsa::twopi(),orsa::twopi());
               const double eclipticLongitude = (tmp_eclipticLongitude>orsa::pi()) ? (tmp_eclipticLongitude-orsa::twopi()) : (tmp_eclipticLongitude);
               const double eclipticLatitude  = theta-theta_sun;
               const double EL = eclipticLongitude;
               const double EB = eclipticLatitude;
            */
            
            if ( (z_H == 160) ||
                 (z_H == 180) ) {
                const double H = z_H*grain_H;
                const double G = 0.15;
                // apparent magnitude
                const double V = apparentMagnitude(H,
                                                   G,
                                                   phaseAngle,
                                                   orb2obs.length(),
                                                   orb2sun.length());
                // detection efficiency
                const double eta = skyCoverage->eta(V,U,AM,GB,GL);
                
                if (!finite(prob_one_more*eta)) {
                    ORSA_DEBUG("skipping: detected non-finite value (prob_one_more: %g   eta: %g   prob_one_more*eta: %g)",prob_one_more,eta,prob_one_more*eta);
                    continue;
                }
                
                xVector[0] = ra*orsa::radToDeg()/15.0;
                xVector[1] = dec*orsa::radToDeg();
                
                if (z_H == 160) {
                    plotStats_sky_NEO_H16->insert(xVector, prob_one_more*eta);
                    // include PHO later...
                    // plotStats_sky_PHO_H16->insert(xVector, eta_PHO);
                }
                
                if (z_H == 180) {
                    plotStats_sky_NEO_H18->insert(xVector, prob_one_more*eta);
                    // include PHO later...
                    // plotStats_sky_PHO_H18->insert(xVector, eta_PHO);
                }

                // now xy (ecliptic plane)
                
                xVector[0] = orsa::FromUnits(orbitPosition_epoch.getX(),orsa::Unit::AU,-1);
                xVector[1] = orsa::FromUnits(orbitPosition_epoch.getY(),orsa::Unit::AU,-1);
                
                if (z_H == 160) {
                    plotStats_xy_NEO_H16->insert(xVector, prob_one_more*eta);
                    // include PHO later...
                    // plotStats_xy_PHO_H16->insert(xVector, eta_PHO);
                }
                
                if (z_H == 180) {
                    plotStats_xy_NEO_H18->insert(xVector, prob_one_more*eta);
                    // include PHO later...
                    // plotStats_xy_PHO_H18->insert(xVector, eta_PHO);
                }
                
                
            }
            
#warning use more general break for NEO and PHO?
            if (sky_N_NEO == target_num) {
                // some output, then break iteration
                break;
            }
        }
        if (sky_N_NEO != target_num) {
            ORSA_DEBUG("left before reaching %i objects",target_num);
        }
        
        {
            // frequent, intermediate output
            if ((callID%100)==0) {
                // ORSA_DEBUG("callID: %Zi -- writing output...",callID.get_mpz_t());
                char filename[1024];
                // 
                sprintf(filename,"tmp_inspect_sky_NEO_H16.dat");
                writeOutputFile(filename, plotStats_sky_NEO_H16, var_ra, var_dec, 1, true);   
                //
                sprintf(filename,"tmp_inspect_sky_NEO_H18.dat");
                writeOutputFile(filename, plotStats_sky_NEO_H18, var_ra, var_dec, 1, true);
                // 
                sprintf(filename,"tmp_inspect_xy_NEO_H16.dat");
                writeOutputFile(filename, plotStats_xy_NEO_H16, var_xy, var_xy, 1);   
                //
                sprintf(filename,"tmp_inspect_xy_NEO_H18.dat");
                writeOutputFile(filename, plotStats_xy_NEO_H18, var_xy, var_xy, 1);
            }
        }
        
    }
    
    ++callID;
    
    return 0;
}

int main(int argc, char ** argv) {
    
    /* 
       {
       // testing only
       unsigned int N_low, N_mean, N_high;
       const double      CL = 0.95;
       const unsigned int R = 15;
       const double       p = 0.75;
       const bool     cache = true; // use cache table in factorial (may use too much memory)
       if (!BinomialConfidenceInterval(N_low, N_mean, N_high, CL, R, p, cache)) {
       ORSA_DEBUG("problems with BinomialConfidenceInterval...");
       }
       exit(0);
       }
    */
    
    if (argc != 3) {
        ORSA_DEBUG("Usage: %s <sqlite-merged-db> <JD>",argv[0]);
        exit(0);
    }
    
    orsa::Debug::instance()->initTimer();
    
    ORSA_DEBUG("process ID: %i",getpid());
    
    // epoch for sky projection
    // global var
    epoch = orsaSolarSystem::FromTimeScale(orsaSolarSystem::julianToTime(atof(argv[2])),
                                           orsaSolarSystem::TS_UTC);
    ORSA_DEBUG("inspect at epoch [below]");
    orsaSolarSystem::print(epoch);
    
    orsaSPICE::SPICE::instance()->loadKernel("de405.bsp");
    
    bg = new orsa::BodyGroup;
    
    // SUN
    sun   = SPICEBody("SUN",orsaSolarSystem::Data::MSun());
    bg->addBody(sun.get());
    
    // EARTH
    earth = SPICEBody("EARTH",orsaSolarSystem::Data::MEarth());
    bg->addBody(earth.get());
    
    // MOON
    moon  = SPICEBody("MOON",orsaSolarSystem::Data::MMoon());
    bg->addBody(moon.get());
    
    earthOrbit.compute(earth.get(),sun.get(),bg.get(),epoch);
    
#warning change random seed
    rnd = new orsa::RNG(2352351);
    
    {
        FitFileDataElement e;
        if (!readFitFile(e,"fit.dat")) {
            ORSA_DEBUG("cannot open file [fit.dat]");
            exit(0); 
        }
        
        skyCoverage = new SkyCoverage;
        
        skyCoverage->obscode  = "500"; // 500=Geocentric
        skyCoverage->epoch    = epoch;
        //
        skyCoverage->V_limit  = e.V_limit;
        skyCoverage->eta0_V   = e.eta0_V;
        skyCoverage->c_V      = e.c_V;
        skyCoverage->w_V      = e.w_V;
        //
        skyCoverage->U_limit  = e.U_limit;
        skyCoverage->w_U      = e.w_U;
        //
        skyCoverage->peak_AM  = e.peak_AM;
        skyCoverage->scale_AM = e.scale_AM;
        skyCoverage->shape_AM = e.shape_AM;
        //
        skyCoverage->drop_GB   = e.drop_GB;
        skyCoverage->scale_GB  = e.scale_GB;
        //
        skyCoverage->scale_GL = e.scale_GL;
        skyCoverage->shape_GL = e.shape_GL;
        //
        skyCoverage->V0       = e.V0;
    }
    
    // needed to work with SQLite database
    sqlite3     * db;
    char        * zErr;
    int           rc;
    std::string   sql;
    
    {
        // open database
        // rc = sqlite3_open(argv[1],&db);
        rc = sqlite3_open_v2(argv[1],&db,SQLITE_OPEN_READONLY,NULL);
        //
        if (rc) {
            fprintf(stderr,"Can't open db: %s\n",sqlite3_errmsg(db));
            sqlite3_close(db);
            exit(0);
        }
    }
    
#warning remember to enable integrity check again
    if (0) {
        // integrity_check
        
        int nrows, ncols;
        char * * sql_result;
        do {
            rc = sqlite3_get_table(db,"pragma integrity_check",&sql_result,&nrows,&ncols,&zErr);
            if (rc==SQLITE_BUSY) {
                ORSA_DEBUG("database busy, retrying...");
                usleep(100000);
            }
        } while (rc==SQLITE_BUSY);
        if (rc != SQLITE_OK) {
            if (zErr != NULL) {
                ORSA_DEBUG("SQL error: %s\n",zErr); 
                sqlite3_free(zErr);
                sqlite3_close(db);
                exit(0);
            }
        }
        //
        bool integrity_check_passed=false;
        if (nrows==1) {
            if (ncols==1) {
                if (sql_result[1]==std::string("ok")) {
                    integrity_check_passed=true;
                }
            }
        }
        if (!integrity_check_passed) {
            ORSA_DEBUG("SQLite problem: integrity_check failed\n"); 
            sqlite3_close(db);
            exit(0); 
        }
        sqlite3_free_table(sql_result);
    }
    
    {
        // now the real work
        
        // vars
        const double a_min  = 1.90;
        const double a_max  = 2.30;
        const double a_step = 0.05;
        osg::ref_ptr<LinearVar> var_a = new LinearVar(a_min,a_max,a_step);
        //
        const double e_min  = 0.00;
        const double e_max  = 1.00;
        const double e_step = 0.05;
        osg::ref_ptr<LinearVar> var_e = new LinearVar(e_min,e_max,e_step);
        //
        const double i_min  =  0.00;
        const double i_max  = 90.00;
        const double i_step =  5.00;
        osg::ref_ptr<LinearVar> var_i = new LinearVar(i_min,i_max,i_step);
        //
        const double L_min  =   0.00;
        const double L_max  = 360.00;
        const double L_step =  30.00;
        osg::ref_ptr<LinearVar> var_L = new LinearVar(L_min,L_max,L_step);
        
        // sky ra,dec [global]
        var_ra  = new LinearVar(  0.0,24.0,0.25);
        var_dec = new LinearVar(-90.0,90.0,5.0);
        
        // xy, in AU
        var_xy  = new LinearVar(-10.0, 10.0, 0.1);
        
        // a,e
        std::vector< osg::ref_ptr<Var> > varDefinition_ae;
        varDefinition_ae.push_back(var_a.get());
        varDefinition_ae.push_back(var_e.get());
        //
        plotStats_ae_NEO_H16 = new PlotStats_WS(varDefinition_ae);
        plotStats_ae_PHO_H16 = new PlotStats_WS(varDefinition_ae);
        //
        plotStats_ae_NEO_H18 = new PlotStats_WS(varDefinition_ae);
        plotStats_ae_PHO_H18 = new PlotStats_WS(varDefinition_ae);
        
        // a,i
        std::vector< osg::ref_ptr<Var> > varDefinition_ai;
        varDefinition_ai.push_back(var_a.get());
        varDefinition_ai.push_back(var_i.get());
        //
        plotStats_ai_NEO_H16 = new PlotStats_WS(varDefinition_ai);
        plotStats_ai_PHO_H16 = new PlotStats_WS(varDefinition_ai);
        //
        plotStats_ai_NEO_H18 = new PlotStats_WS(varDefinition_ai);
        plotStats_ai_PHO_H18 = new PlotStats_WS(varDefinition_ai);
        
        // a,L
        std::vector< osg::ref_ptr<Var> > varDefinition_aL;
        varDefinition_aL.push_back(var_a.get());
        varDefinition_aL.push_back(var_L.get());
        //
        plotStats_aL_NEO_H16 = new PlotStats_WS(varDefinition_aL);
        plotStats_aL_PHO_H16 = new PlotStats_WS(varDefinition_aL);
        //
        plotStats_aL_NEO_H18 = new PlotStats_WS(varDefinition_aL);
        plotStats_aL_PHO_H18 = new PlotStats_WS(varDefinition_aL);
        
        // sky
        std::vector< osg::ref_ptr<Var> > varDefinition_sky;
        varDefinition_sky.push_back(var_ra.get());
        varDefinition_sky.push_back(var_dec.get());
        //
        plotStats_sky_NEO_H16 = new PlotStats_Sum(varDefinition_sky);
        plotStats_sky_PHO_H16 = new PlotStats_Sum(varDefinition_sky);
        //
        plotStats_sky_NEO_H18 = new PlotStats_Sum(varDefinition_sky);
        plotStats_sky_PHO_H18 = new PlotStats_Sum(varDefinition_sky);

        // xy (ecliptic plane)
        std::vector< osg::ref_ptr<Var> > varDefinition_xy;
        varDefinition_xy.push_back(var_xy.get());
        varDefinition_xy.push_back(var_xy.get());
        //
        plotStats_xy_NEO_H16 = new PlotStats_Sum(varDefinition_xy);
        plotStats_xy_PHO_H16 = new PlotStats_Sum(varDefinition_xy);
        //
        plotStats_xy_NEO_H18 = new PlotStats_Sum(varDefinition_xy);
        plotStats_xy_PHO_H18 = new PlotStats_Sum(varDefinition_xy);
        
        {
            char sql_line[1024];
            sprintf(sql_line,"SELECT * FROM grid");
            do {
                rc = sqlite3_exec(db,sql_line,(&inspectCallback),NULL,NULL);
                if (rc==SQLITE_BUSY) {
                    ORSA_DEBUG("database busy, retrying...");
                    usleep(100000);
                }
            } while (rc==SQLITE_BUSY);
            if (rc != SQLITE_OK) {
                if (zErr != NULL) {
                    ORSA_DEBUG("SQL error: %s\n",zErr);
                    sqlite3_free(zErr);
                    sqlite3_close(db);
                    exit(0);
                }
            }
        }
        
        // numbers for Nsub
#warning these can change at each run
        // no sub_a since it is always plotted
        const unsigned int sub_e = 20;
        const unsigned int sub_i = 18;
        const unsigned int sub_node = 12;
        const unsigned int sub_peri = 12;
        const unsigned int sub_M    = 12;
        
        // normalize plotStats_sky_* counts by the area of each bin
        // #warning THIS IS NEEDED!!!
        // normalization moved into writeoutput... functions
        
        // time to write output files
        char filename[1024];
        
        // H16
        sprintf(filename,"%s_inspect_ae_NEO_H16.dat",argv[1]);
        writeOutputFile(filename, plotStats_ae_NEO_H16, var_a, var_e, sub_i*sub_node*sub_peri*sub_M);
        //
        sprintf(filename,"%s_inspect_ae_PHO_H16.dat",argv[1]);
        writeOutputFile(filename, plotStats_ae_PHO_H16, var_a, var_e, sub_i*sub_node*sub_peri*sub_M);
        //
        sprintf(filename,"%s_inspect_ai_NEO_H16.dat",argv[1]);
        writeOutputFile(filename, plotStats_ai_NEO_H16, var_a, var_i, sub_e*sub_node*sub_peri*sub_M);
        //
        sprintf(filename,"%s_inspect_ai_PHO_H16.dat",argv[1]);
        writeOutputFile(filename, plotStats_ai_PHO_H16, var_a, var_i, sub_e*sub_node*sub_peri*sub_M);
        //
        sprintf(filename,"%s_inspect_aL_NEO_H16.dat",argv[1]);
        writeOutputFile(filename, plotStats_aL_NEO_H16, var_a, var_L, sub_e*sub_i*sub_node*sub_peri*3);
        //
        sprintf(filename,"%s_inspect_aL_PHO_H16.dat",argv[1]);
        writeOutputFile(filename, plotStats_aL_PHO_H16, var_a, var_L, sub_e*sub_i*sub_node*sub_peri*3);
        //
        sprintf(filename,"%s_inspect_sky_NEO_H16.dat",argv[1]);
        writeOutputFile(filename, plotStats_sky_NEO_H16, var_ra, var_dec, 1, true);
        //
        sprintf(filename,"%s_inspect_sky_PHO_H16.dat",argv[1]);
        writeOutputFile(filename, plotStats_sky_PHO_H16, var_ra, var_dec, 1, true);
        //
        sprintf(filename,"%s_inspect_xy_NEO_H16.dat",argv[1]);
        writeOutputFile(filename, plotStats_xy_NEO_H16, var_ra, var_dec, 1);
        //
        sprintf(filename,"%s_inspect_xy_PHO_H16.dat",argv[1]);
        writeOutputFile(filename, plotStats_xy_PHO_H16, var_ra, var_dec, 1);
        
        // H18        
        sprintf(filename,"%s_inspect_ae_NEO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_ae_NEO_H18, var_a, var_e, sub_i*sub_node*sub_peri*sub_M);
        //
        sprintf(filename,"%s_inspect_ae_PHO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_ae_PHO_H18, var_a, var_e, sub_i*sub_node*sub_peri*sub_M);
        //
        sprintf(filename,"%s_inspect_ai_NEO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_ai_NEO_H18, var_a, var_i, sub_e*sub_node*sub_peri*sub_M);
        //
        sprintf(filename,"%s_inspect_ai_PHO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_ai_PHO_H18, var_a, var_i, sub_e*sub_node*sub_peri*sub_M);
        //
        sprintf(filename,"%s_inspect_aL_NEO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_aL_NEO_H18, var_a, var_L, sub_e*sub_i*sub_node*sub_peri*3);
        //
        sprintf(filename,"%s_inspect_aL_PHO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_aL_PHO_H18, var_a, var_L, sub_e*sub_i*sub_node*sub_peri*3);
        //
        sprintf(filename,"%s_inspect_sky_NEO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_sky_NEO_H18, var_ra, var_dec, 1, true);
        //
        sprintf(filename,"%s_inspect_sky_PHO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_sky_PHO_H18, var_ra, var_dec, 1, true);
        //
        sprintf(filename,"%s_inspect_xy_NEO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_xy_NEO_H18, var_ra, var_dec, 1);
        //
        sprintf(filename,"%s_inspect_xy_PHO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_xy_PHO_H18, var_ra, var_dec, 1);
    }
    
    sqlite3_close(db);
    
    return 0;
}

