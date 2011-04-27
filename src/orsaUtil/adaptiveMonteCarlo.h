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
        AdaptiveMonteCarlo() : osg::Referenced() {
            
        }
    protected:
        virtual ~AdaptiveMonteCarlo() { }
    public:
        // two different integers between 0 and N-1, with z1 < z2
        /* static bool z1z2(size_t & z1, size_t & z2, const size_t & N) {
           if (N<2) return false;
           // don't use ___z1 and ___z2 in code below, but z1 and z2 
           size_t ___z1 = orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform_int(N);
           size_t ___z2;
           do { ___z2 = orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform_int(N); }
           while (___z1 == ___z2);
           // save them as z1 < z2 to ensure that vectors...[z1] are set before vectors...[z2]
           z1 = std::min(___z1,___z2);
           z2 = std::max(___z1,___z2);
           return true;
           }
        */
        // for this one, we don't force z1 < z2, so it can be z1 > z2
        /* static bool z1z2(size_t & z1, size_t & z2, const size_t & N) {
           if (N<2) return false;
           z1 = orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform_int(N);
           do { z2 = orsa::GlobalRNG::instance()->rng()->gsl_rng_uniform_int(N); }
           while (z1 == z2);
           return true;
           }
        */
        /* 
           virtual void sample(ElementType & e1,
           ElementType & e2,
           const T & intervalVector,
           const size_t & z1,
           const size_t & z2) const {
           }
        */
    public:
        // most implementations will want to rewrite sample(...)
        virtual bool sample(ElementTypeVector & ev) const {
            for (size_t k=0; k<N; ++k) {
                ev[k].position = intervalVector[k]->sample();
#warning maybe add an updateData() to update the auxiliary data...
                intervalVector[k]->updateLevel(ev[k]);
            }
            return true;
        }
        
    public:
        virtual bool run(const T & iv,
                         const size_t & maxIter) {
            intervalVector=iv;
            N=intervalVector.size();
            N.lock();
            std::vector<size_t> intervalVectorOldSize; intervalVectorOldSize.resize(N); for (size_t k=0; k<N; ++k) { intervalVectorOldSize[k] = 0; }
            doAbort=false;
            size_t iter=0;
            // size_t ___z1=0, ___z2=0; // init once to make compiler happy
            while (iter<maxIter) {
                ++iter;
                if (doAbort) break;
                
                // don't use ___z1 and ___z2 in code below, but z1 and z2
                /* if (!z1z2(___z1,___z2,N)) { ORSA_DEBUG("problems..."); }
                   const size_t z1 = ___z1; const size_t z2 = ___z2;
                   // ORSA_DEBUG("z1: %i   z2: %i",z1,z2);
                   */
                
                /* ElementType e1, e2;
                   e1.position = 0.0;
                   e1.level = 0.0;
                   // sample(e1, e2, intervalVector, z1, z2);
                   */
                
                ElementTypeVector ev;
                ev.resize(N);
                sample(ev);

                // insert
                for (size_t k=0; k<N; ++k) {
                    intervalVector[k]->insert(ev[k]);
                }
                
                // every so often, check if a better nominal solution has been found ... old minRMS
                
                
                // every so often, try to decrease the threshold level
#warning should test this only if there has been an insert...       
                // another test: can decreast threshold level?
                if (iter%100==0) {
                    for (size_t k=0; k<N; ++k) {
#warning should test if the adaptive interval actually shrinks when reducing level and at the same time reducing the number of points...
                        if (intervalVector[k]->size()==intervalVectorOldSize[k]) continue;
                        intervalVectorOldSize[k]=intervalVector[k]->size();
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
                            ORSA_DEBUG("reducing threshold for range %3i: %8.2f -> %8.2f   size: %3i -> %3i   range: [%7.3f;%7.3f]",
                                       k,
                                       intervalVector[k]->getThreshold(),
                                       testThreshold,
                                       intervalVector[k]->size(),
                                       testSize,
                                       (*testSampleMin),
                                       (*testSampleMax));
                            intervalVector[k]->updateThresholdLevel(testThreshold);
                        }
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
            return true;
        }
    protected:
        virtual void singleStepDone() const { }
    public:
        virtual void abort() const { doAbort=true; }
    private:
        mutable bool doAbort;
    public:
        T intervalVector;
    public:
        orsa::Cache<size_t> N;
    protected:
        // osg::ref_ptr<orsa::RNG> rng;
    };
    
}; // namespace orsaUtil

#endif // _ORSA_UTIL_ADAPTIVE_MONTECARLOL_H_
