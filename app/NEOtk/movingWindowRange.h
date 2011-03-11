#ifndef _MOVING_WINDOW_RANGE_H_
#define _MOVING_WINDOW_RANGE_H_

#include <deque>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/util.h>

class MovingWindowRange : public osg::Referenced {
public:
    class MovingWindowRangeElement {
    public:
        MovingWindowRangeElement() { }
    public:
        double value; // sampled range value
        double level; // chi-squared 
    public:
        void print() const {
            ORSA_DEBUG("value: %g  level: %g",value,level);
        }
    public:
        bool operator < (const MovingWindowRangeElement & rhs) const {
            return (level < rhs.level);
        }
    };
public:
    /* bool MovingWindowRangeElementCompare(MovingWindowRangeElement * e1,
       MovingWindowRangeElement * e2) {
       return (e1->level < e2->level);
       }
    */
public:
    MovingWindowRange(const double & initialMin,
                      const double & initialMax,
                      const double & samplingTollerance,
                      const unsigned int numPoints,
                      const double & targetLevel,
                      const   bool & forcePositiveRange,
                      const    int & randomSeed) :
        osg::Referenced(),
        min(std::min(initialMin,initialMax)),
        max(std::max(initialMin,initialMax)),
        eps(samplingTollerance),
        maxSize(numPoints),
        target(targetLevel),
        forcePositive(forcePositiveRange),
        rnd(new orsa::RNG(randomSeed))
        {
            if ( (forcePositive && min<0.0) ||
                 (eps <= 0.0) ||
                 (maxSize == 0) ) {
                ORSA_DEBUG("problems");
                exit(0);
            }
        }
protected:
    ~MovingWindowRange() { }
    
public:
    double sample() const {
        return ( (forcePositive?std::max(0.0,min-(max-min)*eps):min-(max-min)*eps) +
                 (max-min)*(1+2*eps)*rnd->gsl_rng_uniform() );
    }
public:
    void insert(const MovingWindowRangeElement & e) {
        
#warning TEST ONLY
       return;
        
        /* ORSA_DEBUG("called insert of:");
           e.print();
        */
        //
        bool doInsert=false;
        if (movingLevel.isSet()) {
            if (e.level < movingLevel.getRef()) {
                doInsert=true;
            }
        } else {
            doInsert=true;
        }
        if (doInsert) {
            data.push_back(e);
            update();
        }
    }
protected:
    void update() {
        // remove point(s) with highest level, if above target level
        // sort(data.begin(),data.end(),MovingWindowRangeElementCompare);
        sort(data.begin(),data.end());
        while (data.size() > maxSize) {
            if (data[data.size()-1].level > target) {
                data.pop_back();
            } else {
                break;
            } 
        }
        movingLevel = data[data.size()-1].level;
        // update min and max
        if (data.size() == maxSize) {
            min = max = data[0].value;
            for (unsigned int k=1; k<data.size(); ++k) {
                if (data[k].value < min) min = data[k].value;
                if (data[k].value > max) max = data[k].value;
            }
        }
    }
public:
    void print() const {
        // for debugging purposes only
        if (movingLevel.isSet()) {
            ORSA_DEBUG("min: %g  max: %g  level: %g",
                       FromUnits(min,orsa::Unit::AU,-1),
                       FromUnits(max,orsa::Unit::AU,-1),
                       movingLevel.getRef());
        } else {
            ORSA_DEBUG("min: %g  max: %g",
                       FromUnits(min,orsa::Unit::AU,-1),
                       FromUnits(max,orsa::Unit::AU,-1));
        }
        /* for (unsigned int k=0; k<data.size(); ++k) {
           ORSA_DEBUG("data[%02i] level: %g  value: %g",k,data[k].level, FromUnits(data[k].value,orsa::Unit::AU,-1));
           }
        */
    }
public:
    typedef std::deque<MovingWindowRangeElement> DataType;
protected:
    DataType data;
public:
    const DataType & getData() const { return data; }
public:
    double getMin() const { return min; }
    double getMax() const { return max; }
protected:
    double min, max;
    const double eps;
    const unsigned int maxSize;
    const double target;
    const bool forcePositive;
protected:
    orsa::Cache<double> movingLevel;
protected:
    osg::ref_ptr<orsa::RNG> rnd;
};

#endif // _MOVING_WINDOW_RANGE_H_
