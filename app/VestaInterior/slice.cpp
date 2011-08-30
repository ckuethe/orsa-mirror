#include <orsa/massDistribution.h>
#include <orsa/matrix.h>
#include <orsa/paulMoment.h>
#include <orsa/util.h>
#include <orsa/vector.h>
#include <orsa/multifit.h>
#include <orsa/chebyshev.h> 
#include <orsa/debug.h>
#include <orsa/print.h>
#include <orsa/statistic.h>

#include <orsaPDS/RadioScienceGravity.h>
#include "CubicChebyshevMassDistribution.h"
#include "simplex.h"

#include "vesta.h"
#include "gaskell.h"

#include <qd/dd_real.h>
#include <qd/qd_real.h>

#include "dislin.h"

using namespace std;
using namespace orsa;

/*** CHOOSE ONE ***/
// typedef double simplex_T;
// typedef mpf_class simplex_T;
typedef dd_real simplex_T;
// typedef qd_real simplex_T;

#warning how to write this using the typedef inside the class?
template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;

/***/

template <typename T> class BinStats : public osg::Referenced {  
    // first, the classes to handle linear and logarithmic variables
public:
    class Var : public osg::Referenced {
    public:
        Var() : osg::Referenced() { }
    protected:
        ~Var() { }
    public:
        virtual size_t size() const = 0;
        virtual size_t bin(const double x) const = 0;
        virtual double binStart(const size_t bin) const = 0;
        virtual double binStop(const size_t bin) const = 0;
        double binCenter(const size_t bin) const {
            return (0.5*(binStart(bin)+binStop(bin)));
        }
    };
public:	
    class LinearVar : public Var {
    public:
        LinearVar(const double startValue,
                  const double stopValue,
                  const double incrValue) : 
            Var(),
            start(startValue),
            stop(stopValue),
            incr(incrValue) {
        }
    public:
        size_t size() const {
            return (size_t)ceil((stop-start)/incr);
        }
        size_t bin(const double x) const {
            const size_t retVal = (size_t)((x-start)/incr);
            // ORSA_DEBUG("x: %g\tbin: %i\tstart: %g\tincr: %g",x,retVal,start,incr);
            // the correct check follows; think twice before changing it!
            if (retVal>=size()) return ((size_t)-1);
            return retVal;
        }
        double binStart(const size_t bin) const {
            return (start+bin*incr);
        }
        double binStop(const size_t bin) const {
            return (start+(bin+1)*incr);
        }
    public:
        const double start, stop, incr;
    };
public:	
    class LogarithmicVar : public Var {
    public:
        LogarithmicVar(const double startValue,
                       const double stopValue,
                       const double factorValue) : 
            Var(),
            start(startValue),
            stop(stopValue),
            factor(factorValue) {
        }
    public:
        size_t size() const {
            return (size_t)ceil(log(stop/start)/log(factor));
        }
        size_t bin(const double x) const {
            const size_t retVal = (size_t)(log(x/start)/log(factor));
            // the correct check follows; think twice before changing it!
            if (retVal>=size()) return ((size_t)-1);
            return retVal;
        }
        double binStart(const size_t bin) const {
            return (start*orsa::int_pow(factor,bin));
        }
        double binStop(const size_t bin) const {
            return (start*orsa::int_pow(factor,bin+1));
        }
    public:
        const double start, stop, factor;
    };
    
public:
    BinStats(const std::vector< osg::ref_ptr<Var> > & varDefinition) :
        osg::Referenced(),
        var(varDefinition) {
        totalSize=1;
        for (unsigned int k=0; k<varDefinition.size(); ++k) {
            totalSize *= varDefinition[k]->size();
        }  
        // ORSA_DEBUG("totalSize: %Zi",totalSize.get_mpz_t());
    }
protected:
    ~BinStats() { }    
public:
    bool bin(std::vector<size_t> & binVector,
             const std::vector<double> & xVector) const {
        if (xVector.size() != var.size()) {
            ORSA_DEBUG("dimension mismatch");
            return false;
        }
        binVector.resize(var.size());
        for (unsigned int k=0; k<var.size(); ++k) {
            binVector[k] = var[k]->bin(xVector[k]);
            if (binVector[k] == ((size_t)-1)) {
                // out of boundaries
                return false;
            }	
        }
        return true;
    }
public:
// from binVector to index
    mpz_class index(const std::vector<size_t> & binVector) const {
        mpz_class idx=0;
        mpz_class mul=1;
        for (unsigned int k=0; k<var.size(); ++k) {
            idx += binVector[k] * mul;
            mul *= var[k]->size();
        }
        // ORSA_DEBUG("index: %i",idx);
        return idx;
    }
public:
// from index to binVector
    std::vector<size_t> bin(mpz_class index) const {
        mpz_class mul=totalSize;
        std::vector<size_t> binVector;
        binVector.resize(var.size());
        {
            unsigned int k=var.size();
            do {
                --k;
                mul /= var[k]->size();
                binVector[k] = mpz_class(index/mul).get_si();
                index -= binVector[k]*mul;
            } while (k>0);
        }
        return binVector;
    }
public:
// vector of the center of each bin
    std::vector<double> binCenterVector(const mpz_class index) const {
        std::vector<size_t> binVector = bin(index);
        std::vector<double> xVector;
        xVector.resize(var.size());
        for (unsigned int k=0; k<var.size(); ++k) {
            xVector[k] = var[k]->binCenter(binVector[k]);
        }
        return xVector;
    }
public:
    mpz_class size() const { return totalSize; }
public:
    const T * stats(const mpz_class index) {
        typename DataType::iterator it = data.find(index);
        if (it != data.end()) {
            return (*it).second.get();
        } else {
            return 0;
        }
    }
protected:
    const std::vector< osg::ref_ptr<Var> > & var;
public:
    typedef typename std::map< mpz_class, osg::ref_ptr<T> > DataType;
protected:
    DataType data;
public:
    const DataType & getData() const {
        return data;
    }
protected:
    mpz_class totalSize;
};

typedef orsa::Statistic<double> PlotStatsElement;

class PlotStats : public BinStats<PlotStatsElement> {
public:
    PlotStats(const std::vector< osg::ref_ptr<Var> > & varDefinition) :
        BinStats<PlotStatsElement>(varDefinition) { }
public:
    bool insert(const std::vector<double> & xVector,
                const double & val) {
        if (xVector.size() != var.size()) {
            ORSA_DEBUG("dimension mismatch");
            return false;
        }
        std::vector<size_t> binVector;
        if (!bin(binVector,xVector)) {
            return false;
        }   
        const mpz_class idx = index(binVector);
        if (data[idx].get()==0) {
            // lazy allocation
            data[idx] = new PlotStatsElement;
        }
        data[idx]->insert(val);

        /* if (xVector.size() > 0) ORSA_DEBUG("xVector[0] = %g",xVector[0]);
           if (xVector.size() > 1) ORSA_DEBUG("xVector[1] = %g",xVector[1]);
           ORSA_DEBUG("insert val: %g",val);
        */
        
        return true;         
    }
};


/***/

int main(int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    // QD
    unsigned int oldcw;
    fpu_fix_start(&oldcw);
    
    //ORSA_DEBUG("current mpf precision: %i",mpf_get_default_prec());
    mpf_set_default_prec(128);
    // mpf_set_default_prec(256);
    // mpf_set_default_prec(512);
    // ORSA_DEBUG("updated mpf precision: %i",mpf_get_default_prec());
    
    if (argc != 5) {
        printf("Usage: %s <plate-model-file> <R0_km> <CCMDF-file> <limit-delta-density_g_cm3>\n",argv[0]);
        exit(0);
    }
    
    const std::string plateModelFile = argv[1];
    const double plateModelR0 = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const std::string CCMDF_filename = argv[3];
    const double limitDeltaDensity = orsa::FromUnits(orsa::FromUnits(atof(argv[4]),orsa::Unit::GRAM),orsa::Unit::CM,-3);
    
    if ( (plateModelR0 <= 0.0) ||
         (limitDeltaDensity <= 0.0) ){
        ORSA_DEBUG("invalid input...");
        exit(0);
    }
    
    // safer over NFS
    sqlite3_vfs_register(sqlite3_vfs_find("unix-dotfile"), 1);
    
    osg::ref_ptr<GaskellPlateModel> shapeModel = new GaskellPlateModel;
    if (!shapeModel->read(plateModelFile)) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    CubicChebyshevMassDistributionFile::DataContainer CCMDF;
    CubicChebyshevMassDistributionFile::read(CCMDF,CCMDF_filename,limitDeltaDensity);
    
    
    /*** plotting ***/
    
    // plot range in km
    const double x_step =    5.00;
    const double x_min  = -300.00;
    const double x_max  =  300.00;
    //
    const double y_step =    5.00;
    const double y_min  = -300.00;
    const double y_max  =  300.00;
    
    std::vector< osg::ref_ptr<PlotStats::Var> > varDefinition;
    //
    osg::ref_ptr<PlotStats::LinearVar> var_x = new PlotStats::LinearVar(x_min,x_max+2*x_step,x_step);
    varDefinition.push_back(var_x.get());
    //
    osg::ref_ptr<PlotStats::LinearVar> var_y = new PlotStats::LinearVar(y_min,y_max+2*y_step,y_step);
    varDefinition.push_back(var_y.get());
    
    osg::ref_ptr<PlotStats> plotStats = 
        new PlotStats(varDefinition);
    
    {
        // random points on given plane
        std::deque<orsa::Vector> rv;
        while (rv.size() < 100000) {
            const orsa::Vector v(orsa::FromUnits(x_min+(x_max-x_min)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform(),orsa::Unit::KM),
                                 orsa::FromUnits(0,orsa::Unit::KM),
                                 orsa::FromUnits(y_min+(y_max-y_min)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform(),orsa::Unit::KM));
            if (shapeModel->isInside(v)) {
                rv.push_back(v);
            }
        }
        std::vector<double> xVector;
        xVector.resize(varDefinition.size());
        CubicChebyshevMassDistributionFile::DataContainer::const_iterator it_CCMDF = CCMDF.begin();
        while (it_CCMDF != CCMDF.end()) {
            osg::ref_ptr<CubicChebyshevMassDistribution> md =
                new CubicChebyshevMassDistribution((*it_CCMDF).coeff,
                                                   (*it_CCMDF).densityScale,
                                                   (*it_CCMDF).R0);
            std::deque<orsa::Vector>::const_iterator it_rv = rv.begin();
            while (it_rv != rv.end()) {
#warning double-check why adding step here...
                xVector[0] = orsa::FromUnits((*it_rv).getX(),orsa::Unit::KM,-1) + x_step;
                xVector[1] = orsa::FromUnits((*it_rv).getZ(),orsa::Unit::KM,-1) + y_step;
                plotStats->insert(xVector,
                                  orsa::FromUnits(orsa::FromUnits(md->density((*it_rv)),orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
                ++it_rv;
            }
            ++it_CCMDF;
        }
    }
    
    const double empty_mesh_val=-1000;
    float * mesh;
    const size_t meshSize = plotStats->size().get_si();
    //
    {
        mesh = (float*)calloc(meshSize,sizeof(float));
        std::vector<double> xVector;
        xVector.resize(varDefinition.size());
        for (unsigned j=0; j<var_x->size(); ++j) {
            for (unsigned k=0; k<var_y->size(); ++k) {
                const unsigned int mesh_id = j*var_y->size()+k;
                xVector[0] = x_min+x_step*(j+0.5);
                xVector[1] = y_min+y_step*(k+0.5);
                std::vector<size_t> binVector;
                if (plotStats->bin(binVector,xVector)) {
                    const PlotStatsElement * e =  plotStats->stats(plotStats->index(binVector));
                    if (e) {
#warning choose what to plot here...
                        if (e->average() > 0.0) {
                            // mesh[mesh_id] = pow10(e->average());
                            mesh[mesh_id] = e->average();
                            // mesh[mesh_id] = log10(e->average());
                        } else {
                            // mesh[mesh_id] = 1e-20;
                            // mesh[mesh_id] = -1000;
                            mesh[mesh_id] = empty_mesh_val;
                         }
                    } else {
                        mesh[mesh_id] = empty_mesh_val;
                    }
                } else {
                    mesh[mesh_id] = empty_mesh_val;
                }
            }
        }
        
        // OLD
        /* const PlotStats::DataType & plotData = plotStats->getData();
           PlotStats::DataType::const_iterator it = plotData.begin();
           while (it != plotData.end()) {
           mesh[(*it).first.get_si()] = (*it).second->average();
           // ORSA_DEBUG("index: %i   val: %g",(*it).first.get_si(),(*it).second->average());
           ++it;
           }
        */
    }
    
    float mesh_min =  1.0e100;
    float mesh_max = -1.0e100;
    for (unsigned int k=0; k<meshSize; ++k) {
        // ORSA_DEBUG("mesh[%06i] = %g",k,mesh[k]);
        if (mesh[k]==empty_mesh_val) continue;
        if (mesh[k]<mesh_min) mesh_min=mesh[k];
        if (mesh[k]>mesh_max) mesh_max=mesh[k];
    }
    
    
    
    /*** DISLIN ***/
    page(4000,4000);
    pagmod("LAND");
    
    // output file name
    /* char plotFilename[1024];
       sprintf(plotFilename,"%s.fit.pdf",basename.c_str());
       setfil(plotFilename);
    */
    // new files overwrite old ones
    filmod("DELETE");
    
    // metafl("POST");
    metafl("PSCL");
    // metafl("PNG");
    // metafl("XWIN"); x11mod("STORE"); clrmod("FULL");
    // metafl("CONS");
    // metafl("TIFF");
    
    pagmod("LAND");
    // page(3600,1400);
    // pagfll(255);
    // winsiz(800,600);
    
    // background color
    scrmod("REVERS");
    
    disini();
    // pagera();
    hwfont();
    simplx();
    // triplx();
    // helve();
    // helves();
    // winfnt();
    // disalf();
    // psfont("AvantGarde-Book");
    // color("fore");
    
    penwid(1.0);
    height(50); // text height
    
    // paghdr("","",2,0);
    
    /* axspos(200,1300);
       axslen(1350,1050);
    */
    axspos( 500,3500);
    axslen(2500,2500);
    
    // select a color table
    setvlt("RAIN"); // TEMP,GREY,RGREY,VGA,RAIN,SPEC...
    // setvlt("TEMP"); // TEMP,GREY,RGREY,VGA,RAIN,SPEC...
    // setvlt("RGREY"); // TEMP,GREY,RGREY,VGA,RAIN,SPEC...
    // setvlt("RAIN");
    
    hwmode("ON","LINE");
    
    texmod("ON"); // TeX text
    
    // NEOs
    // titlin("NEOs In-Field Probability for 703",4);
    // titlin("NEOs In-Field Probability for G96",4);
    // titlin("H=18 NEOs Detection Efficiency for 703",4);
    // titlin("H=18 NEOs Detection Efficiency for G96",4);
    // titlin("H=18 NEOs Observation Probability for 703",4);
    // titlin("H=18 NEOs Observation Probability for G96",4);
    // 
    // PHOs
    // titlin("PHOs In-Field Probability for 703",4);
    // titlin("PHOs In-Field Probability for G96",4);
    // titlin("H=18 PHOs Detection Efficiency for 703",4);
    // titlin("H=18 PHOs Detection Efficiency for G96",4);
    // titlin("H=18 PHOs Observation Probability for 703",4);
    // titlin("H=18 PHOs Observation Probability for G96",4);
    //
    // misc
    // titlin("H=20 PHOs Observation Probability for 703",4);
    // titlin("H=20 PHOs Observation Probability for G96",4);
    // titlin("H=22 PHOs Observation Probability for 703",4);
    // titlin("H=22 PHOs Observation Probability for G96",4);
    
    // titlin("3-D Colour Plot of the Function",2);
    // titlin("F(X,Y) = 2 * SIN(X) * SIN(Y)",4);
    // titlin("Saturn Trojans Predictor",4);

    // titlin("TITLE HERE",4);
    
    // name("initial semi-major axis [AU]","x");
    // name("libration amplitude","x");
    // name("eccentricity","y");
    // name("Z-axis","z");

    // a,e
    name("[km]","x");
    name("[km]","y");
    name("density [g/cm$^3$]","z");
    // a,i
    /* name("Semi-Major Axis [AU]","x");
       name("Inclination [deg]","y");
    */
    // a,L
    /* name("Semi-Major Axis [AU]","x");
       name("True Longitude [deg]","y");
    */
    //
    // name("Probability","z");
    // name("Detection Efficiency","z");
    // name("Probability","z");
    
    // name("Long.","x");
    // name("Lat.","y");
    
    // axsscl("log","z");
    // labels("float","y");  

    setgrf("NAME","NAME","NONE","NONE");
    
    intax();

    // a,e
    /* autres(var_x->size()+3,
       var_y->size()+3);
    */
    /* autres(var_x->size()+1,
       var_y->size()+2);
    */
    autres(var_x->size(),
           var_y->size());
    // a,i
    /* autres(var_x->size()+3,
       var_y->size()+3);
    */
    // a,L
    /* autres(var_x->size()+3,
       var_y->size()+1);
    */
    
    // labtyp("vert","z"); // vertical labels for z axis

    // a,e
    //
    // NEOs
    // const double z_min=1e-4; const double z_max=1e-1;
    // const double z_min=1e-3; const double z_max=1e-1;
    // const double z_min=1e-5; const double z_max=1e-3;
    //
    // PHOs
    // const double z_min=1e-4; const double z_max=1e-1;
    // const double z_min=1e-3; const double z_max=1e0;
    // const double z_min=1e-6; const double z_max=1e-3;
    //
    // other... LINEAR 0->1
    // const double z_min=1.0e-5; const double z_max=1.0e-3;
    const double z_min=2.0; const double z_max=5.0;
    //
    // a,i
    // const double z_min=1e-5; const double z_max=1e-3;
    //
    // a,L
    // const double z_min=1e-5; const double z_max=1e-3;
    //
    // 
    {
        // bound z
        for (unsigned int k=0; k<meshSize; ++k) {
            if (mesh[k]==empty_mesh_val) continue;
            if (mesh[k]<z_min) mesh[k]=z_min;
            if (mesh[k]>z_max) mesh[k]=z_max;
        }
    }
    // a,e
    /* {
       axsscl("log","Z");
       labels("log","Z");
       digits(1,"X"); 
       digits(1,"Y");
       digits(0,"Z");
       ticks(1,"X");
       ticks(1,"Y");
       ticks(5,"Z");
       frame(5); // frame thickness
       graf3(x_min-x_step/2,x_max+x_step/2,x_min,25.0,
       y_min-y_step/2,y_max+y_step/2,y_min,25.0,
       log10(z_min),log10(z_max),log10(z_min),1.0);
       crvmat(mesh,var_x->size(),var_y->size(),1,1);
       }
    */
    // a,e
    {
        digits(0,"X"); 
        digits(0,"Y");
        digits(0,"Z");
        ticks(2,"X");
        ticks(2,"Y");
        ticks(2,"Z");
        // axsscl("log","z");
        // labels("log","z");
        frame(5); // frame thickness
        graf3(x_min-0.5*x_step,x_max+0.5*x_step,x_min,100.0,
              y_min-0.5*y_step,y_max+0.5*y_step,y_min,100.0,
              z_min,z_max,z_min,1.0);
        crvmat(mesh,var_x->size(),var_y->size(),1,1);
        if (0) {
            // q=1.3 line
            const size_t nray=1000;
            float xray[nray],yray[nray];
            for (size_t k=0; k<nray; ++k) {
                xray[k] = 1.9+((2.3-1.9)*k)/nray;
                yray[k] = 1.0-(1.3/xray[k]);
                ORSA_DEBUG("ray-> %g %g",xray[k],yray[k]);
            }
            graf(x_min-0.5*x_step,x_max+0.5*x_step,0,0,
                 y_min-0.5*y_step,y_max+0.5*y_step,0,0);
            lintyp(2);
            curve(xray,yray,nray);
        }
    }
    // a,i
    /* {
       digits(1,"X"); 
       digits(0,"Y");
       digits(0,"Z");
       ticks(1,"X");
       ticks(1,"Y");
       ticks(5,"Z");
       axsscl("log","z");
       labels("log","z");
       frame(5); // frame thickness
       graf3(x_min-x_step/2,x_max+x_step/2,1.0,0.2,
       y_min-y_step/2,y_max+y_step/2,0.0,15.0,
       log10(z_min),log10(z_max),log10(z_min),1);
       crvmat(mesh,var_x->size(),var_y->size(),1,1);
       }
    */
    // a,L
    /* {
       digits(1,"X"); 
       digits(0,"Y");
       digits(0,"Z");
       ticks(1,"X");
       ticks(3,"Y");
       ticks(5,"Z");
       axsscl("log","z");
       labels("log","z");
       frame(5); // frame thickness
       graf3(x_min-x_step/2,x_max+x_step/2,1.0,0.2,
       y_min-y_step/2,y_max+y_step/2,y_min,90.0,
       log10(z_min),log10(z_max),log10(z_min),1);
       crvmat(mesh,var_x->size(),var_y->size(),1,1);
       }
    */
    
    
    /* 
       if (1) {
       // contour plot, level curves
       digits(1,"contour");
       labels("FLOAT","CONTUR");
       double T=eta_min+round_eta;
       while (T<eta_max) {
       conmat((float *)mesh,NX,NY,T);
       T += round_eta;
       }
       }
    */
    
    // title only
    vkytit(-50); // title closer to plot
    height(50); // text height
    title();
    
    disfin();
    
    return 0;
}
