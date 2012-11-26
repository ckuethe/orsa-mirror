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
        printf("Usage: %s <plate-model-file> <R0_km> <CCMDF-file> <XY|XZ|YZ>\n",argv[0]);
        exit(0);
    }
    
    const std::string plateModelFile = argv[1];
    const double plateModelR0 = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const std::string CCMDF_filename = argv[3];
    const std::string mode_str = argv[4];
    
    enum MODE {XY,XZ,YZ};
    MODE mode;
    if (mode_str == "XY") {
        mode = XY;
    } else if (mode_str == "XZ") {
        mode = XZ;
    } else if (mode_str == "YZ") {
        mode = YZ;
    } else {
        ORSA_DEBUG("mode [%s] not recognized",mode_str.c_str());
        exit(0);
    }
    
    if (plateModelR0 <= 0.0){
        ORSA_DEBUG("invalid input...");
        exit(0);
    }
    
    // safer over NFS
    // sqlite3_vfs_register(sqlite3_vfs_find("unix-dotfile"), 1);
    
    osg::ref_ptr<GaskellPlateModel> shapeModel = new GaskellPlateModel;
    if (!shapeModel->read(plateModelFile)) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    
    CubicChebyshevMassDistributionFile::DataContainer CCMDF;
    CubicChebyshevMassDistributionFile::read(CCMDF,CCMDF_filename);
    if (CCMDF.size() == 0) exit(0);
    
    /*** plotting ***/

    // plot range in km
    double xs,xm,xM;
    double ys,ym,yM;
    switch (mode) {
        case XY:
            xs=0.30; xm=-80.0; xM=100.0;
            ys=0.30; ym=-70.0; yM= 70.0;
            break;
        case XZ:
            xs=0.30; xm=-80.0; xM=100.0;
            ys=0.30; ym=-70.0; yM= 70.0;
            break;
        case YZ:
            xs=0.30; xm=-70.0; xM= 70.0;
            ys=0.30; ym=-70.0; yM= 70.0;
            break;
    }
    const double x_step = xs;
    const double x_min  = xm;
    const double x_max  = xM;
    //
    const double y_step = ys;
    const double y_min  = ym;
    const double y_max  = yM;
    
    // plot range in km
    /* const double x_step =    5.00;
       const double x_min  = -300.00;
       const double x_max  =  300.00;
       //
       const double y_step =    5.00;
       const double y_min  = -300.00;
       const double y_max  =  300.00;
    */
    //
    /* const double x_step =    0.30;
       const double x_min  = - 80.00;
       const double x_max  =  100.00;
       //
       const double y_step =    0.30;
       const double y_min  = - 50.00;
       const double y_max  =   50.00;
    */
    
    std::vector< osg::ref_ptr<PlotStats::Var> > varDefinition;
    //
    osg::ref_ptr<PlotStats::LinearVar> var_x = new PlotStats::LinearVar(x_min,x_max+2*x_step,x_step);
    varDefinition.push_back(var_x.get());
    //
    osg::ref_ptr<PlotStats::LinearVar> var_y = new PlotStats::LinearVar(y_min,y_max+2*y_step,y_step);
    varDefinition.push_back(var_y.get());
    
    /* osg::ref_ptr<PlotStats> plotStats = 
       new PlotStats(varDefinition);
    */
    
    FILE * fp_xyz = fopen("slice.xyz","w");
    
    {
        ORSA_DEBUG("sampling...");
        // random points on given plane
        std::deque<orsa::Vector> rv_in;
        std::deque<orsa::Vector> rv_out;
        {
            double x = x_min + 0.5*x_step;
            while (x<x_max) {
                double y = y_min + 0.5*y_step;
                while (y<y_max) {

                    double vx,vy,vz;
                    switch (mode) {
                        case XY:
                            vx=x; vy=y; vz=0;
                            break;
                        case XZ:
                            vx=x; vy=0; vz=y;
                            break;
                        case YZ:
                            vx=0; vy=x; vz=y;
                            break;
                    }
                    const orsa::Vector v(orsa::FromUnits(vx,orsa::Unit::KM),
                                         orsa::FromUnits(vy,orsa::Unit::KM),
                                         orsa::FromUnits(vz,orsa::Unit::KM));
                    
                    // XZ
                    /* const orsa::Vector v(orsa::FromUnits(x,orsa::Unit::KM),
                       orsa::FromUnits(0,orsa::Unit::KM),
                       orsa::FromUnits(y,orsa::Unit::KM));
                    */
                    
                    // XY
                    /* const orsa::Vector v(orsa::FromUnits(x,orsa::Unit::KM),
                       orsa::FromUnits(y,orsa::Unit::KM),
                       orsa::FromUnits(0,orsa::Unit::KM));
                    */
                    
                    if (shapeModel->isInside(v)) {
                        rv_in.push_back(v);
                    } else {
                        rv_out.push_back(v);
                    }
                    
                    y += y_step;
                }
                x += x_step;
            }
        }
        //
        ORSA_DEBUG("sampling... done.");
        //
        std::vector<double> xVector;
        xVector.resize(varDefinition.size());
        // double toKM = orsa::FromUnits(1,orsa::Unit::KM,-1);
        CubicChebyshevMassDistributionFile::DataContainer::const_iterator it_CCMDF = CCMDF.begin();
        
#warning EDIT: PLOT LAST ONE ONLY...
        it_CCMDF = CCMDF.end();
        --it_CCMDF;
        
        while (it_CCMDF != CCMDF.end()) {
            osg::ref_ptr<CubicChebyshevMassDistribution> md = CCMD(*it_CCMDF);
            std::deque<orsa::Vector>::const_iterator it_rv_in = rv_in.begin();
            while (it_rv_in != rv_in.end()) {
                // #warning double-check why adding step here...
                
                switch (mode) {
                    case XY:
                        xVector[0] = orsa::FromUnits((*it_rv_in).getX(),orsa::Unit::KM,-1);
                        xVector[1] = orsa::FromUnits((*it_rv_in).getY(),orsa::Unit::KM,-1);
                        break;
                    case XZ:
                        xVector[0] = orsa::FromUnits((*it_rv_in).getX(),orsa::Unit::KM,-1);
                        xVector[1] = orsa::FromUnits((*it_rv_in).getZ(),orsa::Unit::KM,-1);
                        break;
                    case YZ:
                        xVector[0] = orsa::FromUnits((*it_rv_in).getY(),orsa::Unit::KM,-1);
                        xVector[1] = orsa::FromUnits((*it_rv_in).getZ(),orsa::Unit::KM,-1);
                        break;
                }
                
                
                // XZ
                /* xVector[0] = orsa::FromUnits((*it_rv_in).getX(),orsa::Unit::KM,-1);
                   xVector[1] = orsa::FromUnits((*it_rv_in).getZ(),orsa::Unit::KM,-1);
                */
                
                // XY
                /* xVector[0] = orsa::FromUnits((*it_rv_in).getX(),orsa::Unit::KM,-1);
                   xVector[1] = orsa::FromUnits((*it_rv_in).getY(),orsa::Unit::KM,-1);
                */
                
                fprintf(fp_xyz,"%g %g %g\n",xVector[0],xVector[1],orsa::FromUnits(orsa::FromUnits(md->density((*it_rv_in)),orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
                
                /* plotStats->insert(xVector,
                   orsa::FromUnits(orsa::FromUnits(md->density((*it_rv_in)),orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
                */
                
                ++it_rv_in;
            }

            // append zeros outside body
            std::deque<orsa::Vector>::const_iterator it_rv_out = rv_out.begin();
            while (it_rv_out != rv_out.end()) {
                // #warning double-check why adding step here...

                switch (mode) {
                    case XY:
                        xVector[0] = orsa::FromUnits((*it_rv_out).getX(),orsa::Unit::KM,-1);
                        xVector[1] = orsa::FromUnits((*it_rv_out).getY(),orsa::Unit::KM,-1);
                        break;
                    case XZ:
                        xVector[0] = orsa::FromUnits((*it_rv_out).getX(),orsa::Unit::KM,-1);
                        xVector[1] = orsa::FromUnits((*it_rv_out).getZ(),orsa::Unit::KM,-1);
                        break;
                    case YZ:
                        xVector[0] = orsa::FromUnits((*it_rv_out).getY(),orsa::Unit::KM,-1);
                        xVector[1] = orsa::FromUnits((*it_rv_out).getZ(),orsa::Unit::KM,-1);
                        break;
                }
                
                
                // XZ
                /* xVector[0] = orsa::FromUnits((*it_rv_out).getX(),orsa::Unit::KM,-1);
                   xVector[1] = orsa::FromUnits((*it_rv_out).getZ(),orsa::Unit::KM,-1);
                */
                
                // XY
                /* xVector[0] = orsa::FromUnits((*it_rv_out).getX(),orsa::Unit::KM,-1);
                   xVector[1] = orsa::FromUnits((*it_rv_out).getY(),orsa::Unit::KM,-1);
                */
                
                fprintf(fp_xyz,"%g %g %g\n",xVector[0],xVector[1],0.0);
                
                ++it_rv_out;
            }

            
            ++it_CCMDF;
        }
    }
    
    fclose(fp_xyz);
    
    ORSA_DEBUG("done.");
    
    return 0;
}
