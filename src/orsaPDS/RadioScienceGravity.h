#ifndef _ORSA_PDS_RADIO_SCIENCE_GRAVITY_H_
#define _ORSA_PDS_RADIO_SCIENCE_GRAVITY_H_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <string>

namespace orsaPDS {
    
    class RadioScienceGravityData : public osg::Referenced {
    public:
        double R0, GM, sigmaGM;
        int degree, order;
        int normalizationState;
        int numberOfCoefficients;
        double referenceLongitude, referenceLatitude;
    };
    
    class RadioScienceGravityFile : public osg::Referenced {
    public:
        // read RECORD_BYTES and FILE_RECORDS from the label of the PDS file
        RadioScienceGravityFile(const std::string & fileName,
                                const size_t & RECORD_BYTES_,
                                const size_t & FILE_RECORDS_);
    public:
        const RadioScienceGravityData * getData() const { return data.get(); }
    protected:
        osg::ref_ptr<RadioScienceGravityData> data;
    protected:
        bool readD(double & d) const;
        bool readI(int & i) const;
        bool readS(std::string & s) const;
    protected:
        void skipToNextRow() const;
        
    protected:
        const size_t RECORD_BYTES;
        const size_t FILE_RECORDS;
        
    protected:
        FILE * fp;
    };
    
} // namespace orsaPDS

#endif // _ORSA_PDS_RADIO_SCIENCE_GRAVITY_H_
