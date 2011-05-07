#ifndef _INSIDE_VESTA_H_
#define _INSIDE_VESTA_H_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/util.h>
#include <orsa/vector.h>
#include <orsa/massDistribution.h>

// Each parameter of the model
class ModelParameter : public osg::Referenced {
public:
    ModelParameter(const double & min_,
                   const double & max_) :
        osg::Referenced(),
        min(std::min(min_,max_)),
        range(fabs(max_-min_)) { }
protected:
    virtual ~ModelParameter() { }
public:
    double sample() const {
        return min+range*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform();
    }
protected:
    const double min, range;
};

typedef ModelParameter Par;

// All the parameters in the model in the interior structure of Vesta
// keep all vars in sync!
class Model {
public:
    // global
    osg::ref_ptr<Par> totalMass;
    osg::ref_ptr<Par> totalVolume;
    // core
    osg::ref_ptr<Par> coreDensity;
    osg::ref_ptr<Par> coreCenterX;
    osg::ref_ptr<Par> coreCenterY;
    osg::ref_ptr<Par> coreCenterZ;
    osg::ref_ptr<Par> coreRadiusX; // X-Y-Z dimension of core BEFORE ROTATION
    osg::ref_ptr<Par> coreRadiusY;
    osg::ref_ptr<Par> coreRadiusZ;
    // osg::ref_ptr<Par> 
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
        _coreVolume((4.0/3.0)*orsa::pi()*_val.coreRadiusX*_val.coreRadiusY*_val.coreRadiusZ),
        _mantleDensity((_val.totalMass-_coreVolume*_val.coreDensity)/(_val.totalVolume-_coreVolume)) {
        if (_coreVolume > 0.0) {
            _am2 = orsa::int_pow((*_val.coreRadiusX),-2);
            _bm2 = orsa::int_pow((*_val.coreRadiusY),-2);
            _cm2 = orsa::int_pow((*_val.coreRadiusZ),-2);
        }
    }
protected:
    ~ModelMassDistribution() { }
public:    
    double density(const orsa::Vector & v) const {
        if (_coreVolume > 0.0) {
            if ( (orsa::square(v.getX()-_val.coreCenterX)*_am2 +
                  orsa::square(v.getY()-_val.coreCenterY)*_bm2 +
                  orsa::square(v.getZ()-_val.coreCenterZ)*_cm2) <= 1.0) {
                // inside core
                return _val.coreDensity;
            } else {
                return _mantleDensity;
            }
        } else {
            return _mantleDensity;
        }
    }
public:
    const Model::Values _val;
    const double _coreVolume;
    const double _mantleDensity;
protected:
    orsa::Cache<double> _am2, _bm2, _cm2;
};

#endif // _INSIDE_VESTA_H_
