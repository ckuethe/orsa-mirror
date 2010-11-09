#ifndef _INSIDE_VESTA_H_
#define _INSIDE_VESTA_H_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/util.h>
#include <orsa/vector.h>
#include <orsa/massDistribution.h>

// A singletone class for the random number generator
class GlobalRNG {   
public:
    static GlobalRNG * instance() {
        if (_instance == 0) {
            _instance = new GlobalRNG;
        }
        return _instance;
    }
protected:
    GlobalRNG() {
        const pid_t pid = getpid();
        ORSA_DEBUG("RNG seed: %i",pid);
        rng = new orsa::RNG(pid);
    }
public:
    virtual ~GlobalRNG() {
        _instance = 0;
    }
protected:
    static GlobalRNG * _instance;
protected:
    osg::ref_ptr<orsa::RNG> rng;
public:
    double gsl_rng_uniform() const {
        return rng->gsl_rng_uniform();
    }
};

// Each parameter of the model
class ModelParameter : public osg::Referenced {
public:
    ModelParameter(const double & min_,
                   const double & max_) :
        osg::Referenced(),
        min(min_),
        range(max_-min_) { }
protected:
    virtual ~ModelParameter() { }
public:
    double sample() const {
        return min+range*GlobalRNG::instance()->gsl_rng_uniform();
    }
protected:
    const double min, range;
};

typedef ModelParameter Par;

// All the parameters in the model in the interior structure of Vesta
// keep all vars in sync!
class Model {
public:
    osg::ref_ptr<Par> totalMass;
    osg::ref_ptr<Par> totalVolume;
    osg::ref_ptr<Par> coreDensity;
    osg::ref_ptr<Par> coreCenterX;
    osg::ref_ptr<Par> coreCenterY;
    osg::ref_ptr<Par> coreCenterZ;
    osg::ref_ptr<Par> coreRadiusX; // X-Y-Z dimension of core BEFORE ROTATION
    osg::ref_ptr<Par> coreRadiusY; 
    osg::ref_ptr<Par> coreRadiusZ;
    // no core rotation for now
    // core-mantle interface
    // no mantle density, that's a dependant variable, to conserve total mass
    
public:
    class Values {
    public:
        orsa::Cache<double> totalMass;
        orsa::Cache<double> totalVolume;
        orsa::Cache<double> coreDensity;
        orsa::Cache<double> coreCenterX;
        orsa::Cache<double> coreCenterY;
        orsa::Cache<double> coreCenterZ;
        orsa::Cache<double> coreRadiusX;
        orsa::Cache<double> coreRadiusY; 
        orsa::Cache<double> coreRadiusZ;
    };
public:
    Values sample() const {
        Values val;
        val.totalMass   = totalMass->sample();
        val.totalVolume = totalVolume->sample();
        val.coreDensity = coreDensity->sample();
        val.coreCenterX = coreCenterX->sample();
        val.coreCenterY = coreCenterY->sample();
        val.coreCenterZ = coreCenterZ->sample();
        val.coreRadiusX = coreRadiusX->sample();
        val.coreRadiusY = coreRadiusY->sample(); 
        val.coreRadiusZ = coreRadiusZ->sample();
        return val;
    }
};

class ModelMassDistribution : public orsa::MassDistribution {
public:
    ModelMassDistribution(const Model::Values & val) :
        orsa::MassDistribution(),
        _val(val),
        _am2(orsa::int_pow(_val.coreRadiusX.getRef(),-2)),
        _bm2(orsa::int_pow(_val.coreRadiusY.getRef(),-2)),
        _cm2(orsa::int_pow(_val.coreRadiusZ.getRef(),-2)),
        _coreVolume((4.0/3.0)*orsa::pi()*_val.coreRadiusX.getRef()*_val.coreRadiusY.getRef()*_val.coreRadiusZ.getRef()),
        _mantleDensity((_val.totalMass.getRef()-_coreVolume*_val.coreDensity.getRef())/(_val.totalVolume.getRef()-_coreVolume)) { }
protected:
    ~ModelMassDistribution() { }
public:    
    double density(const orsa::Vector & v) const {
        if ( (orsa::square(v.getX()-_val.coreCenterX.getRef())*_am2 +
              orsa::square(v.getY()-_val.coreCenterY.getRef())*_bm2 +
              orsa::square(v.getZ()-_val.coreCenterZ.getRef())*_cm2) <= 1.0) {
            // inside core
            return _val.coreDensity.getRef();
        } else {
            return _mantleDensity;
        }
    }
public:
    const Model::Values _val;
    const double _am2, _bm2, _cm2;
    const double _coreVolume;
    const double _mantleDensity;
};

#endif // _INSIDE_VESTA_H_
