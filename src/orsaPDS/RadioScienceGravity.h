#ifndef _ORSA_PDS_RADIO_SCIENCE_GRAVITY_H_
#define _ORSA_PDS_RADIO_SCIENCE_GRAVITY_H_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <string>

namespace orsaPDS {
    
    class RadioScienceGravityData : public osg::Referenced {
        
    };
    
    class RadioScienceGravityFile : public osg::Referenced {
    public:
        RadioScienceGravityFile(const std::string & fileName);
    public:
        const RadioScienceGravityData * getData() const { return data.get(); }
    protected:
        osg::ref_ptr<RadioScienceGravityData> data;
    };
    
} // namespace orsaPDS

#endif // _ORSA_PDS_RADIO_SCIENCE_GRAVITY_H_
