#ifndef _SPLIT_RANGE_H_
#define _SPLIT_RANGE_H_

#include <deque>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/util.h>

// all linear for now, then will split logarithmically

// linear is easier also because we do not need to check if range min and max are positive

class RangeElement {
public:
    double min, max;
    double weight;
};

// assume each RangeElement covers an unique range...

class Range : public osg::Referenced {
public:
    Range(const double & min,
          const double & max,
          const unsigned int initialSplit,
          const int & randomSeed) : rnd(new orsa::RNG(randomSeed)) {
        for (unsigned int k=0; k<initialSplit; ++k) {
            RangeElement re;
            re.min = min+(max-min)*k/initialSplit;
            re.max = min+(max-min)*(k+1)/initialSplit;
            re.weight = 1.0;
            data.push_back(re);
        }
        normalizeWeights();
    }
protected:
    virtual ~Range() { }
protected:
    void normalizeWeights() {
        /* const double sum=sumWeights();
           for (unsigned int k=0; k<data.size(); ++k) {
           data[k].weight /= sum;
           }
        */
    }
protected:
    double sumWeights() const {
        double sum=0.0;
        for (unsigned int k=0; k<data.size(); ++k) {
            if (data[k].weight < 0.0) {
                ORSA_DEBUG("problem: negative weight");
            }
            sum += data[k].weight;
        }
        return sum;
    }
public:
    double sample() const {
        const double sum_w = sumWeights();
        const double x = sum_w*rnd->gsl_rng_uniform();
        double runningWeight=0.0;
        for (unsigned int k=0; k<data.size(); ++k) {
            if ( (x>=runningWeight) &&
                 (x<(runningWeight+data[k].weight)) ) {
                return (data[k].min + (data[k].max-data[k].min)*rnd->gsl_rng_uniform());
            }
            runningWeight += data[k].weight;
        }
        ORSA_DEBUG("problems");
        return 0.0;
    }
public:
    void feedback(const double & r,
                  const double & weightFactor) {

#warning TEST ONLY
        return;
        
        for (unsigned int k=0; k<data.size(); ++k) {
            if ( (r>=data[k].min) &&
                 (r< data[k].max) ) {
                data[k].weight *= weightFactor;
#warning set here min and max weight for each split range
                if (data[k].weight <  0.1) data[k].weight =  0.1;
                if (data[k].weight > 10.0) data[k].weight = 10.0;
#warning leak in nearby if reaching maximum?
                if (data[k].weight == 10.0) {
                    if ( (k>0) &&
                         (data[k-1].weight < 1.0) ) {
                        data[k-1].weight = 1.0;
                    }
                    if ( ((k+1)<data.size()) &&
                         (data[k+1].weight < 1.0) ) {
                        data[k+1].weight = 1.0;
                    }
                }
                break;
            }
        }
        normalizeWeights();
    }
public:
    void feedbackAll(const double & weightFactor) {
        
#warning TEST ONLY
        return;
        
        for (unsigned int k=0; k<data.size(); ++k) {
            data[k].weight *= weightFactor;
#warning set here min and max weight for each split range
            if (data[k].weight <  0.1) data[k].weight =  0.1;
            if (data[k].weight > 10.0) data[k].weight = 10.0;
        }
    }
    
public:
    void split(const double & weightThreshold) {
        ORSA_DEBUG("code needed here, for 1 delete and 2 inserts, for each split range...");
    }
public:
    size_t size() const { return data.size(); }
public:
    typedef std::deque<RangeElement> DataType;
protected:
    DataType data;
public:
    const DataType & getData() const { return data; }
    
protected:
    osg::ref_ptr<orsa::RNG> rnd;
};

#endif // _SPLIT_RANGE_H_
