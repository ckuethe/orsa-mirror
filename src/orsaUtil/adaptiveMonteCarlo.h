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
        virtual bool run(std::vector< osg::ref_ptr<orsaUtil::AdaptiveInterval<T> > > & intervalVector,
                         const size_t & maxIter) {
            doAbort=false;
            size.lock();
            size_t iter=0;
            while (iter<maxIter) {
                ++iter;
                
            }
            size.unlock();
            return true;
        }
    protected:
        virtual void singleStepDone() const { }
    public:
        virtual void abort() const { doAbort=true; }
    private:
        mutable bool doAbort;
    public:
        orsa::Cache<size_t> size;
    };
    
}; // namespace orsaUtil

#endif // _ORSA_UTIL_ADAPTIVE_MONTECARLOL_H_
