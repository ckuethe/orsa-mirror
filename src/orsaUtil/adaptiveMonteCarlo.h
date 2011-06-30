#ifndef _ORSA_UTIL_ADAPTIVE_MONTECARLOL_H_
#define _ORSA_UTIL_ADAPTIVE_MONTECARLOL_H_

#include <orsa/interval.h>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/util.h>

#include <orsaUtil/adaptiveInterval.h>

namespace orsaUtil {
    
    template <typename T> class AdaptiveMonteCarlo : public osg::Referenced {
    public:
        typedef typename T::value_type::element_type::AdaptiveIntervalElementType ElementType;
        typedef typename std::vector<ElementType> ElementTypeVector;
    public:
        AdaptiveMonteCarlo() : osg::Referenced() { }
    protected:
        virtual ~AdaptiveMonteCarlo() { }
    public:
        // most derived classes will want to reimplement sample(...), in particular to update all the data within each element
        virtual bool sample(ElementTypeVector & ev) const {
            ++iterCount;
            for (size_t k=0; k<N; ++k) {
                ev[k].position = intervalVector[k]->sample();
                ev[k].position.lock();
                intervalVector[k]->updateLevel(ev[k]);
            }
            return true;
        }
        
    public:
        virtual bool run(const T & iv,
                         const size_t & maxIter) {
            intervalVector=iv;
            maxIterCount=maxIter;
            N=intervalVector.size();
            N.lock();
            std::vector<size_t> intervalVectorOldSize; intervalVectorOldSize.resize(N); for (size_t k=0; k<N; ++k) { intervalVectorOldSize[k] = 0; }
            doAbort=false;
            iterCount=0;
            while (iterCount<maxIterCount) {
                // ++iter; // increase iterCount within sample(...), as sample(...) can require multiple iterations to produce a good sample
                if (doAbort) break;
                
                ElementTypeVector ev;
                ev.resize(N);
                if (!sample(ev)) {
                    ORSA_DEBUG("problems...");
                    return false;
                }
                
                // insert
                for (size_t k=0; k<N; ++k) {
                    intervalVector[k]->insert(ev[k]);
                }
                
                if (iterCount%100==0) {
                    
                    // every so often, check if a better nominal solution has been found ... old minRMS
                    periodicUpdate();
                    
                    // every so often, try to decrease the threshold level
#warning should test this only if there has been an insert...       
                    // another test: can decrease threshold level?
                    for (size_t k=0; k<N; ++k) {
#warning should test if the adaptive interval actually shrinks when reducing level and at the same time reducing the number of points...
                        if (intervalVector[k]->size()==intervalVectorOldSize[k]) continue;
                        // intervalVectorOldSize[k]=intervalVector[k]->size(); // update oldVec later...
                        if (intervalVector[k]->getThreshold() == intervalVector[k]->getTargetThreshold()) continue;
                        const double currentSampleRange = intervalVector[k]->maxSample()-intervalVector[k]->minSample();
#warning parameters...
                        // const size_t minSize = 50;
                        const double thresholdIncreaseFactor = 1.2;
                        // if (intervalVector[k]->size()<minSize) continue;
                        double testThreshold = intervalVector[k]->getTargetThreshold();
                        double testSampleRange = currentSampleRange;
                        orsa::Cache<double> testMin, testMax;
                        orsa::Cache<double> testSampleMin, testSampleMax;
                        size_t testSize = 0;
                        while (1) {
                            if (testThreshold>=intervalVector[k]->getThreshold()) break;
                            testSize = 0;
                            {
                                typename T::value_type::element_type::DataType::DataType::const_iterator it =
                                    intervalVector[k]->getData()->getData().begin();
                                while (it != intervalVector[k]->getData()->getData().end()) {
                                    if ((*it).level < testThreshold) {
                                        testMin.setIfSmaller((*it).position);
                                        testMax.setIfLarger((*it).position);
                                        ++testSize;
                                    }
                                    ++it;
                                }
                            }
                            // test if new interval would be smaller than current one
                            if (testMin.isSet() && testMax.isSet()) {
                                // should use an actual AdaptiveInterval for this?
                                
                                // same as 2*delta in
                                // must use same residual probability as actual intervals
#warning use number different from 2 ?  
                                if (testSize >= 2) {
                                    T::value_type::element_type::sampleRange(testSampleMin,
                                                                             testSampleMax,
                                                                             intervalVector[k]->initialMin,
                                                                             intervalVector[k]->initialMax,
                                                                             testMin,
                                                                             testMax,
                                                                             intervalVector[k]->probability,
                                                                             testSize);
                                    testSampleRange = testSampleMax - testSampleMin;
#warning use a factor to compare ranges? if so, correct every occurrence of testSampleRange...
                                    if (testSampleRange < currentSampleRange) {
                                        break; // done
                                    }
                                }
                            }
                            // exception break: have enough points at target threshold
                            if ( (testThreshold==intervalVector[k]->getTargetThreshold()) && 
                                 (testSize>=intervalVector[k]->targetSamples) ) {
                                ORSA_DEBUG("QUICK EXIT for intervalVector[%i] ****************************",k);
                                break;
                            }
                            testThreshold *= thresholdIncreaseFactor;
                        }
                        bool reduceThreshold=false;
                        if ( (testSampleRange < currentSampleRange) && 
                             (testSize>=2) &&
                             (testThreshold<intervalVector[k]->getThreshold()) ) {
                            reduceThreshold=true;
                        }
                        if ( (testThreshold==intervalVector[k]->getTargetThreshold()) &&
                             (testSize>=intervalVector[k]->targetSamples) ) {
                            reduceThreshold=true;
                        }
                        if (reduceThreshold) {
                            // version with range in AU
                            /* ORSA_DEBUG("reducing threshold for range %3i: %8.2f -> %8.2f   size: %3i -> %3i   range: %7.3f -> %7.3f [AU]",
                               k,
                               intervalVector[k]->getThreshold(),
                               testThreshold,
                               intervalVector[k]->size(),
                               testSize,
                               orsa::FromUnits(currentSampleRange,orsa::Unit::AU,-1),
                               orsa::FromUnits(testSampleRange,orsa::Unit::AU,-1));
                            */
                            ORSA_DEBUG("reducing threshold for interval %3i: %8.2f -> %8.2f   size: %3i -> %3i   range: [%.2g;%.2g]   delta: %.2g",
                                       k,
                                       intervalVector[k]->getThreshold(),
                                       testThreshold,
                                       intervalVector[k]->size(),
                                       testSize,
                                       (*testSampleMin),
                                       (*testSampleMax),
                                       testSampleMax-testSampleMin);
                            intervalVector[k]->updateThresholdLevel(testThreshold);
                        }
                    }
                    
                    for (size_t k=0; k<N; ++k) {
                        intervalVectorOldSize[k]=intervalVector[k]->size();
                    }
                    
                }
                
                
                // check if we can stop now
                {
                    bool canStop=true;
                    for (size_t k=0; k<N; ++k) {
                        if ( (intervalVector[k]->size() >= intervalVector[k]->targetSamples) &&
                             (intervalVector[k]->getThreshold() == intervalVector[k]->getTargetThreshold()) ) {
                            // still can stop, so keep canStop=true;
                        } else {
                            canStop=false;
                        }
                    }
                    if (canStop) break;
                }
            }
            N.unlock();
            ORSA_DEBUG("iterCount: %i",iterCount);
            return true;
        }
    protected:
        virtual void periodicUpdate() const { }
    public:
        virtual void abort() const { doAbort=true; }
    private:
        mutable bool doAbort;
    public:
        T intervalVector;
    public:
        orsa::Cache<size_t> N;
    protected:
        mutable size_t iterCount, maxIterCount;
    };
    
}; // namespace orsaUtil

#endif // _ORSA_UTIL_ADAPTIVE_MONTECARLOL_H_
