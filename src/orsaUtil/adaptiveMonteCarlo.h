#ifndef _ORSA_UTIL_ADAPTIVE_MONTECARLOL_H_
#define _ORSA_UTIL_ADAPTIVE_MONTECARLOL_H_

#include <orsa/interval.h>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/util.h>

#include <orsaUtil/adaptiveInterval.h>

namespace orsaUtil {
    
    class AdaptiveMonteCarlo : public osg::Referenced {
    public:
        AdaptiveMonteCarlo();
    protected:
        virtual ~AdaptiveMonteCarlo();
    public:
        virtual bool run(const size_t & maxIter);
    protected:
        virtual void singleStepDone() const;
    public:
        virtual void abort() const { doAbort=true; }
    private:
        mutable bool doAbort;
    public:
        
    };
    
}; // namespace orsaUtil

#endif // _ORSA_UTIL_ADAPTIVE_MONTECARLOL_H_
