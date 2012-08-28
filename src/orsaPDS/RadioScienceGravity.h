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
        // friend class orsaPDS::RadioScienceGravityFile;
    public:
        double R0, GM, sigmaGM;
        unsigned int degree, order;
        unsigned int normalizationState;
        unsigned int numberOfCoefficients;
        double referenceLongitude, referenceLatitude;
    public:
        // maps the coefficient name to its index
        // keys are: GM, C002000, C002001, S002001, C002002, S002002, ...
    public:
        typedef QHash<QString, unsigned int> HashType;
    public:
        HashType hash;
    public:
        static QString keyC(unsigned int l, unsigned int m);
        static QString keyS(unsigned int l, unsigned int m);
    public:
        unsigned int index(const QString & key) const;
        QString key(const unsigned int & index) const;
    public:
        std::vector< orsa::Cache<double> >                coeff;
        std::vector< std::vector< orsa::Cache<double> > > covar; // triangular
        /* public:
           void resize(const unsigned int degree_) {
           degree = degree_;
           numberOfCoefficients = (degree+1)*(degree+1)-3;
           coeff.resize(numberOfCoefficients);
           covar.resize(numberOfCoefficients);
           for (unsigned int k=0; k<numberOfCoefficients; ++k) {
           covar[k].resize(k+1);
           }
           }
        */
    public:
        double getCoeff(const QString &) const;
        double getCovar(const QString &, const QString &) const;
    public:
        void setCoeff(const QString &, const double &);
        void setCovar(const QString &, const QString &, const double &);
    public:
        gsl_vector * getCoefficientVector() const;
        gsl_matrix * getCovarianceMatrix() const;
        gsl_matrix * getInverseCovarianceMatrix() const;
    };
    
    class RadioScienceGravityFile {
        // read RECORD_BYTES and FILE_RECORDS from the label of the PDS file
        
        // file READ interface
    public:
        static bool read(RadioScienceGravityData * data,
                         const std::string & fileName,
                         const size_t & RECORD_BYTES,
                         const size_t & FILE_RECORDS);
    protected:
        static bool readD(double & d, FILE * fp);
        static bool readI(int & i, FILE * fp);
        static bool readU(unsigned int & u, FILE * fp);
        static bool readS(std::string & s, FILE * fp);
    protected:
        static void skipToNextRow(const size_t & RECORD_BYTES, FILE * fp);
        
        // file WRITE interface
    public:
        static bool write(const RadioScienceGravityData * data,
                          const std::string & fileName,
                          const size_t & RECORD_BYTES,
                          const size_t & FILE_RECORDS);
    protected:
        static bool writeD(const double & d, FILE * fp);
        static bool writeI(const int & i, FILE * fp);
        static bool writeU(const unsigned int & u, FILE * fp);
        static bool writeS(const std::string & s, FILE * fp);
    protected:
        static void padToNextRow(const size_t & RECORD_BYTES, FILE * fp);
    };
    
} // namespace orsaPDS

#endif // _ORSA_PDS_RADIO_SCIENCE_GRAVITY_H_
