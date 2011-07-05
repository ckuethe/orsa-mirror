#ifndef _ORSA_PDS_RADIO_SCIENCE_GRAVITY_H_
#define _ORSA_PDS_RADIO_SCIENCE_GRAVITY_H_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/vector.h>

#include <string>
#include <vector>

#include <QHash>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace orsaPDS {

    class RadioScienceGravityFile;
    
    class RadioScienceGravityData : public osg::Referenced {
        friend class orsaPDS::RadioScienceGravityFile;
    public:
        double R0, GM, sigmaGM;
        unsigned int degree, order;
        unsigned int normalizationState;
        unsigned int numberOfCoefficients;
        double referenceLongitude, referenceLatitude;
    protected:
        // maps the coefficient name to its index
        // keys are: GM, C002000, C002001, S002001, C002002, S002002, ...
    public:
        typedef QHash<QString, unsigned int> HashType;
    protected:
        HashType hash;
    public:
        static QString keyC(unsigned int l, unsigned int m);
        static QString keyS(unsigned int l, unsigned int m);
    protected:
        std::vector< orsa::Cache<double> >                coeff;
        std::vector< std::vector< orsa::Cache<double> > > covar; // triangular
    public:
        double getCoeff(const QString &) const;
        double getCovar(const QString &, const QString &) const;
    public:
        gsl_vector * getCoefficientVector() const;
        gsl_matrix * getCovarianceMatrix() const;
        gsl_matrix * getInverseCovarianceMatrix() const;
    };
    
    class RadioScienceGravityFile : public osg::Referenced {
    public:
        // read RECORD_BYTES and FILE_RECORDS from the label of the PDS file
        RadioScienceGravityFile(const std::string & fileName,
                                const size_t & RECORD_BYTES_,
                                const size_t & FILE_RECORDS_);
    public:
        const RadioScienceGravityData * getData() const { return data.get(); }
    public:
        osg::ref_ptr<RadioScienceGravityData> data;
    protected:
        bool readD(double & d) const;
        bool readI(int & i) const;
        bool readU(unsigned int & u) const;
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
