#ifndef _ORSA_UTIL_ADAPTIVE_INTERVAL_H_
#define _ORSA_UTIL_ADAPTIVE_INTERVAL_H_

#include <orsa/interval.h>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/util.h>

namespace orsaUtil {
    
    // template <typename T> class AdaptiveIntervalElement : public T {
    template <typename T> class AdaptiveIntervalElement {
        // position in the interval
        // level such as chi-squared
    public:
        AdaptiveIntervalElement() {
            data = new T;
        }
    public:
        AdaptiveIntervalElement(const AdaptiveIntervalElement & aie) {
            if (aie.position.isSet()) position = aie.position;
            if (aie.level.isSet()) level = aie.level;
            data = aie.data;
        }
    public:
        const AdaptiveIntervalElement & operator = (const AdaptiveIntervalElement & aie) {
            if (aie.position.isSet()) position = aie.position;
            if (aie.level.isSet()) level = aie.level;
            data = aie.data;
            return (*this);
        }
    public:
        /* AdaptiveIntervalElement(const T * t) {
           data = t;
           }
        */
    public:
        orsa::Cache<double> position;
        mutable orsa::Cache<double> level;
    public:
        mutable osg::ref_ptr<T> data;
    public:
        inline bool operator == (const AdaptiveIntervalElement & rhs) const {
            return (position == rhs.position);
        }
        inline bool operator != (const AdaptiveIntervalElement & rhs) const {
            return (position != rhs.position);
        }
        inline bool operator <  (const AdaptiveIntervalElement & rhs) const {
            return (position <  rhs.position);
        }
        inline bool operator >  (const AdaptiveIntervalElement & rhs) const {
            return (position >  rhs.position);
        }
        inline bool operator <= (const AdaptiveIntervalElement & rhs) const {
            return (position <= rhs.position);
        }
        inline bool operator >= (const AdaptiveIntervalElement & rhs) const {
            return (position >= rhs.position);
        }
    };
    
    template <typename T> class AdaptiveInterval : public osg::Referenced {
    public:
        AdaptiveInterval(const double & min,
                         const double & max,
                         const double & residualProbability, // 1-confidenceLevel
                         const double & initialThresholdLevel,
                         const double & targetThresholdLevel) :
            osg::Referenced(),
            initialMin(std::min(min,max)),
            initialMax(std::max(min,max)),
            probability(residualProbability),
            threshold(initialThresholdLevel),
            targetThreshold(targetThresholdLevel) {
            if ( (probability < 0.0) ||
                 (probability > 1.0) ||
                 (threshold < targetThreshold)) {
                ORSA_DEBUG("problems");
                exit(0);
            }
            data = new DataType;
            data->enableDataStoring();
        }
        
    protected:
        ~AdaptiveInterval() { }
        
    public:
        double minSample() const {
            if (data->size()>=2) {
                return sampleMin;
            } else {
                return initialMin;
            }
        }
        double maxSample() const {
            if (data->size()>=2) {
                return sampleMax;
            } else {
                return initialMax;
            }
        }
        double sample() const {
            // #warning minimum number of points to use (at least 2...)
            if (data->size()>=2) {
                /* 
                   #warning cache values, to speed-up this part 
                   // const double    tmpMin = forcePositive ? std::max(0.0,data->min().position) : data->min().position;
                   // forcePositive test is not needed, because of the stronger check below agains the intial min and max
                   const double    tmpMin = data->min().position;
                   const double    tmpMax = data->max().position;
                   const double        mu = 0.5*(tmpMin+tmpMax);
                   const double     delta = (tmpMax-tmpMin)/(2*pow(probability,1.0/data->size()));
                   const double sampleMin = std::max(mu-delta,initialMin);
                   const double sampleMax = std::min(mu+delta,initialMax);
                */
                return (sampleMin+(sampleMax-sampleMin)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform());
            } else {
                return (initialMin+(initialMax-initialMin)*orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform());
            }
        }
    public:
        virtual void updateLevel(const AdaptiveIntervalElement<T> & e) const = 0;
    public:
        virtual void updateThresholdLevel(const double & newThreshold) {
            // ORSA_DEBUG("NEW threshold: %g",newThreshold);
            threshold = newThreshold;
            // ORSA_DEBUG("size BEFORE refresh: %i",data->size());
            typename DataType::DataType::iterator it = data->getData().begin();
            while (it != data->getData().end()) {
                // (*it).level.reset();
                // updateLevel(*it);
                if ((*it).level < threshold) {
                    // ORSA_DEBUG("keeping level: %g",(*(*it).level));
                    ++it;
                } else {
                    // ORSA_DEBUG("rejecting level: %g",(*(*it).level));
                    it = data->getData().erase(it);
                }
            }
            data->update();
            update();
            // ORSA_DEBUG("size AFTER refresh: %i",data->size());
        }
    public:
        double getThreshold() const { return threshold; } 
        double getTargetThreshold() const { return targetThreshold; } 
    public:
        void insert(const AdaptiveIntervalElement<T> & e) {
            if (!e.level.isSet()) {
                updateLevel(e);
            }
            // ORSA_DEBUG("testing level: %g   threshold: %g",(*e.level),threshold);
            if (e.level < threshold) {
                data->insert(e,false,false);
                data->update();
                update();
            }
        }
    public:
        // call refresh() when vecRMS changes
        void refresh() {
            // ORSA_DEBUG("size BEFORE refresh: %i",data->size());
            typename DataType::DataType::iterator it = data->getData().begin();
            while (it != data->getData().end()) {
                (*it).level.reset();
                updateLevel(*it);
                if ((*it).level < threshold) {
                    ++it;
                } else {
                    it = data->getData().erase(it);
                }
            }
            data->update();
            update();
            // ORSA_DEBUG("size AFTER refresh: %i",data->size());
        }
    public:
        bool reset() {
            return data->reset();
        }
    public:
        void update() {
            if (data->size()>=2) {
                const double    tmpMin = data->min().position;
                const double    tmpMax = data->max().position;
                //
                const double        mu = 0.5*(tmpMin+tmpMax);
                const double     delta = (tmpMax-tmpMin)/(2*pow(probability,1.0/data->size()));
                //
                sampleMin = std::max(mu-delta,initialMin);
                sampleMax = std::min(mu+delta,initialMax);
            }
        }
    protected:
        orsa::Cache<double> sampleMin, sampleMax;
    public:
        size_t size() const { return data->size(); }
    public:
        typedef typename orsaUtil::AdaptiveIntervalElement<T> AdaptiveIntervalElementType;
        typedef typename orsa::Interval<AdaptiveIntervalElementType> DataType;
    protected:
        osg::ref_ptr<DataType> data;
    public:
        const DataType * getData() const { return data.get(); }
        
    protected:
        const double initialMin;
        const double initialMax;
        const double probability;
    protected:
        double threshold;
        const double targetThreshold;
    };

}; // namespace orsaUtil

#endif // _ORSA_UTIL_ADAPTIVE_INTERVAL_H_
