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
            for (unsigned int k=0; k<N; ++k) {
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

                // every so often, check if a better nominal solution has been found ... old minRMS
                
                
                // every so often, try to decrease the threshold level

                
                // check if we can stop now
                
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
