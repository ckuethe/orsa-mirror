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
        AdaptiveMonteCarlo() : osg::Referenced() {
            
        }
    protected:
        virtual ~AdaptiveMonteCarlo() { }
    public:
        // two different integers between 0 and N-1, with z1 < z2
        static bool z1z2(size_t & z1, size_t & z2, const size_t & N) {
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
    public:
        virtual bool run(T & intervalVector,
                         const size_t & maxIter) {
            N=intervalVector.size();
            N.lock();
            doAbort=false;
            size_t iter=0;
            unsigned int ___z1=0, ___z2=0; // init once to make compiler happy
            while (iter<maxIter) {
                ++iter;
                if (doAbort) break;
                
                // don't use ___z1 and ___z2 in code below, but z1 and z2
                if (!z1z2(___z1,___z2,N)) { ORSA_DEBUG("problems..."); }
                const size_t z1 = ___z1; const size_t z2 = ___z2;
                ORSA_DEBUG("z1: %i   z2: %i",z1,z2);
                
                
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
        orsa::Cache<size_t> N;
    protected:
        // osg::ref_ptr<orsa::RNG> rng;
    };
    
}; // namespace orsaUtil

#endif // _ORSA_UTIL_ADAPTIVE_MONTECARLOL_H_
