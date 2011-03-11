#ifndef _ADAPTIVE_INTERVAL_H_
#define _ADAPTIVE_INTERVAL_H_

#include <orsa/interval.h>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/util.h>

template <typename T> class AdaptiveIntervalElement {
    // position in the interval
    // level such as chi-squared
public:
    AdaptiveIntervalElement() { }
public:
    orsa::Cache<double> position;
    orsa::Cache<double> level;
public:
    T t;
public:
    inline bool operator == (const AdaptiveIntervalElement & rhs) const {
        return (position.getRef() == rhs.position.getRef());
    }
    inline bool operator != (const AdaptiveIntervalElement & rhs) const {
        return (position.getRef() != rhs.position.getRef());
    }
    inline bool operator <  (const AdaptiveIntervalElement & rhs) const {
        return (position.getRef() <  rhs.position.getRef());
    }
    inline bool operator >  (const AdaptiveIntervalElement & rhs) const {
        return (position.getRef() >  rhs.position.getRef());
    }
    inline bool operator <= (const AdaptiveIntervalElement & rhs) const {
        return (position.getRef() <= rhs.position.getRef());
    }
    inline bool operator >= (const AdaptiveIntervalElement & rhs) const {
        return (position.getRef() >= rhs.position.getRef());
    }
};

template <typename T> class AdaptiveInterval : public osg::Referenced {
public:
    AdaptiveInterval(const double & min,
                     const double & max,
                     const double & confidenceLevel,
                     const double & thresholdLevel,
#warning remove forcePositive? it is not needed if we constantly check to be within initial range...
                     const   bool & forcePositiveRange,
                     const    int & randomSeed) :
        osg::Referenced(),
        initialMin(std::min(min,max)),
        initialMax(std::max(min,max)),
        probability(1.0-confidenceLevel),
        threshold(thresholdLevel),
        forcePositive(forcePositiveRange),
        rnd(new orsa::RNG(randomSeed))
        {
            if ( (forcePositive && initialMin<0.0) ||
                 (probability < 0.0) ||
                 (probability > 1.0) ) {
                ORSA_DEBUG("problems");
                exit(0);
            }
            data = new DataType;
            data->enableDataStoring();
        }
protected:
    ~AdaptiveInterval() { }
    
      
public:
    double sample() const {
        if (data->size()>=2) {
#warning cache values, to speed-up this part 
            // const double    tmpMin = forcePositive ? std::max(0.0,data->min().position.getRef()) : data->min().position.getRef();
            // forcePositive test is not needed, because of the stronger check below agains the intial min and max
            const double    tmpMin = data->min().position.getRef();
            const double    tmpMax = data->max().position.getRef();
            /* ORSA_DEBUG("tmpMin: %g  tmpMax: %g",
               orsa::FromUnits(tmpMin,orsa::Unit::AU,-1),
               orsa::FromUnits(tmpMax,orsa::Unit::AU,-1));
            */
            const double        mu = 0.5*(tmpMin+tmpMax);
            const double     delta = (tmpMax-tmpMin)/(2*pow(probability,1.0/data->size()));
            const double sampleMin = std::max(mu-delta,initialMin);
            const double sampleMax = std::min(mu+delta,initialMax);
            /* ORSA_DEBUG("sampleMin: %g  sampleMax: %g  mu: %g  delta: %g  probability: %g",
               orsa::FromUnits(sampleMin,orsa::Unit::AU,-1),
               orsa::FromUnits(sampleMax,orsa::Unit::AU,-1),
               orsa::FromUnits(mu,orsa::Unit::AU,-1),
               orsa::FromUnits(delta,orsa::Unit::AU,-1),
               probability);
            */
            return (sampleMin+(sampleMax-sampleMin)*rnd->gsl_rng_uniform());
        } else {
            return (initialMin+(initialMax-initialMin)*rnd->gsl_rng_uniform());
        }
    }
public:
    void insert(const AdaptiveIntervalElement<T> & e) {
        if (e.level.getRef() < threshold) {
            /* ORSA_DEBUG("inserting pos: %g  level: %g",
               orsa::FromUnits(e.position.getRef(),orsa::Unit::AU,-1),
               e.level.getRef());
            */
            data->insert(e,false,false);
        }
    }
public:
    typedef typename orsa::Interval< AdaptiveIntervalElement<T> > DataType;
protected:
    osg::ref_ptr<DataType> data;
public:
    const DataType * getData() const { return data.get(); }
    
protected:
    const double initialMin;
    const double initialMax;
    const double probability;
    const double threshold;
    const   bool forcePositive;
protected:
    osg::ref_ptr<orsa::RNG> rnd;
};

#endif // _ADAPTIVE_INTERVAL_H_
