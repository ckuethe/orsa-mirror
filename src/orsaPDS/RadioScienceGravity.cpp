#include <orsaPDS/RadioScienceGravity.h>

#include <iostream>
#include <cstdio>

#include <orsa/print.h>
#include <orsa/util.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

using namespace orsaPDS;

QString RadioScienceGravityData::keyC(unsigned int l, unsigned int m) {
    QString qs;
    qs.sprintf("C%03i%03i",l,m);
    return qs;
}

QString RadioScienceGravityData::keyS(unsigned int l, unsigned int m) {
    QString qs;
    qs.sprintf("S%03i%03i",l,m);
    return qs;
}

unsigned int RadioScienceGravityData::index(const QString & s) const {
    const HashType::const_iterator it = hash.constFind(s);
    if (it != hash.constEnd()) {
        return it.value();
    } else {
        ORSA_DEBUG("problem: no hash found for [%s]",s.toStdString().c_str());
        return 0;
    }
}

double RadioScienceGravityData::getCoeff(const QString & s) const {
    const HashType::const_iterator it = hash.constFind(s);
    if (it != hash.constEnd()) {
        return coeff[it.value()];
    } else {
        ORSA_DEBUG("problem: no hash found for [%s]",s.toStdString().c_str());
        return 0;
    }
}

double RadioScienceGravityData::getCovar(const QString & s1, const QString & s2) const {
    const HashType::const_iterator it1 = hash.constFind(s1);
    const HashType::const_iterator it2 = hash.constFind(s2);
    if ( (it1 != hash.constEnd()) && (it2 != hash.constEnd()) ) {
        return covar[std::max(it1.value(),it2.value())][std::min(it1.value(),it2.value())];
    } else {
        if (it1 == hash.constEnd()) ORSA_DEBUG("problem: no hash found for [%s]",s1.toStdString().c_str());
        if (it2 == hash.constEnd()) ORSA_DEBUG("problem: no hash found for [%s]",s2.toStdString().c_str());
        return 0;
    }
}

gsl_vector * RadioScienceGravityData::getCoefficientVector() const {
    gsl_vector * mu = gsl_vector_alloc(numberOfCoefficients);
    for (unsigned int k=0; k<numberOfCoefficients; ++k) {
        gsl_vector_set(mu,k,coeff[k]);
    }
    return mu;
}

gsl_matrix * RadioScienceGravityData::getCovarianceMatrix() const {
    gsl_matrix * covm = gsl_matrix_alloc(numberOfCoefficients,numberOfCoefficients);
    for (unsigned int l=0; l<numberOfCoefficients; ++l) {
        for (unsigned int m=0; m<numberOfCoefficients; ++m) { 
            if (m<=l) {
                gsl_matrix_set(covm,l,m,covar[l][m]);
            } else {
                gsl_matrix_set(covm,l,m,covar[m][l]);
            }
        }
    }
    return covm;
}

gsl_matrix * RadioScienceGravityData::getInverseCovarianceMatrix() const {
    // first, find eigenvectors and eigenvalues of covariance matrix
    gsl_eigen_symmv_workspace * workspace = gsl_eigen_symmv_alloc(numberOfCoefficients);
    gsl_matrix * evec = gsl_matrix_alloc(numberOfCoefficients,numberOfCoefficients);
    gsl_vector * eval = gsl_vector_alloc(numberOfCoefficients);
    gsl_eigen_symmv(getCovarianceMatrix(),eval,evec,workspace);
    // check if definite positive
    for (size_t k=0;k<numberOfCoefficients;++k) {
        if (gsl_vector_get(eval,k) <= 0.0) {
            ORSA_DEBUG("problems... eval[%03i] = %g",k,gsl_vector_get(eval,k));
            return 0;
        }
    }
    // make a diagonal matrix with eigenvalues    
    /* gsl_matrix * eval_matrix = gsl_matrix_alloc(numberOfCoefficients,numberOfCoefficients);
       for (size_t i=0;i<numberOfCoefficients;++i) {
       for (size_t j=0;j<numberOfCoefficients;++j) {
       if (i==j) {
       gsl_matrix_set(eval_matrix,i,j,gsl_vector_get(eval,i));
       } else {
       gsl_matrix_set(eval_matrix,i,j,0.0);
       }
       }
       }
    */
    // inverse of eval matrix (diagonal with values 1/eigenval)
    gsl_matrix * inverse_eval_matrix = gsl_matrix_alloc(numberOfCoefficients,numberOfCoefficients);
    for (size_t i=0;i<numberOfCoefficients;++i) {
        for (size_t j=0;j<numberOfCoefficients;++j) {
            if (i==j) {
                gsl_matrix_set(inverse_eval_matrix,i,j,1.0/gsl_vector_get(eval,i));
            } else {
                gsl_matrix_set(inverse_eval_matrix,i,j,0.0);
            }
        }
    }
    gsl_matrix * inv_covm = gsl_matrix_alloc(numberOfCoefficients,numberOfCoefficients);
    gsl_matrix * tmp_matrix = gsl_matrix_alloc(numberOfCoefficients,numberOfCoefficients);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,inverse_eval_matrix,evec,0.0,tmp_matrix);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,evec,tmp_matrix,0.0,inv_covm);
    
    // free all gsl allocated vars, except inv_covm
    gsl_eigen_symmv_free(workspace);    
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
    // gsl_matrix_free(eval_matrix);
    gsl_matrix_free(inverse_eval_matrix);
    // skip inv_covm
    gsl_matrix_free(tmp_matrix);
    
    return inv_covm;
}

RadioScienceGravityFile::RadioScienceGravityFile(const std::string & fileName,
                                                 const size_t & RECORD_BYTES_,
                                                 const size_t & FILE_RECORDS_) :
    RECORD_BYTES(RECORD_BYTES_),
    FILE_RECORDS(FILE_RECORDS_) {
    
    data = new RadioScienceGravityData;
    
    fp = fopen(fileName.c_str(),"r");
    if (!fp) {
        ORSA_DEBUG("cannot open file [%s]",fileName.c_str());
        return;
    }
    
    // read SHBDR_HEADER_TABLE
    {
        readD(data->R0);
        data->R0 = orsa::FromUnits(data->R0,orsa::Unit::KM);
        
        readD(data->GM);
        data->GM = orsa::FromUnits(orsa::FromUnits(data->GM,orsa::Unit::KM,3),orsa::Unit::SECOND,-2);
        
        readD(data->sigmaGM);
        data->sigmaGM = orsa::FromUnits(orsa::FromUnits(data->sigmaGM,orsa::Unit::KM,3),orsa::Unit::SECOND,-2);
        
        readU(data->degree);
        
        readU(data->order);
        
        readU(data->normalizationState);

        readU(data->numberOfCoefficients);
        
        readD(data->referenceLongitude);
        data->referenceLongitude *= orsa::degToRad();
        
        readD(data->referenceLatitude);
        data->referenceLatitude *= orsa::degToRad();
        
        skipToNextRow();
    }
    
    // ORSA_DEBUG("n. coeff: %i",data->numberOfCoefficients);
    
    data->coeff.resize(data->numberOfCoefficients);

    data->covar.resize(data->numberOfCoefficients);
    for (unsigned int k=0; k<data->numberOfCoefficients; ++k) {
        data->covar[k].resize(k+1);
    }
    
    // read SHBDR_NAMES_TABLE
    {
        std::string s;
        for (unsigned int k=0; k<data->numberOfCoefficients; ++k) {
            readS(s);
            // ORSA_DEBUG("name: [%s]",s.c_str());
            data->hash[QString(s.c_str())] = k;
        }
        
        skipToNextRow();
    }
    
    // read SHBDR_COEFFICIENTS_TABLE
    {
        double d;
        for (unsigned int k=0; k<data->numberOfCoefficients; ++k) {
            readD(d);
            // ORSA_DEBUG("d: [%g]",d);
            if (k == data->index("GM")) {
                // GM is in km^2/s^2
                d = orsa::FromUnits(orsa::FromUnits(d,orsa::Unit::KM,3),orsa::Unit::SECOND,-2);
            }
            data->coeff[k] = d;
        }
        
        skipToNextRow();
    }
    
    // read SHBDR_COVARIANCE_TABLE
    {
        double d;
        for (unsigned int k=0; k<data->numberOfCoefficients; ++k) {
            for (unsigned int r=0; r<=k; ++r) {
                readD(d);
                // ORSA_DEBUG("d: [%g]",d);
                data->covar[k][r] = d;
            }
        }
        
        skipToNextRow();
    }
    
    fclose(fp);
}


bool RadioScienceGravityFile::readD(double & d) const {
    const bool retVal = (1 == fread(&d,sizeof(double),1,fp));
    if (!retVal) {
        ORSA_DEBUG("problems...");
    }
    return retVal;
}

bool RadioScienceGravityFile::readI(int & i) const {
    const bool retVal = (1 == fread(&i,sizeof(int),1,fp));
    if (!retVal) {
        ORSA_DEBUG("problems...");
    }
    return retVal;
}

bool RadioScienceGravityFile::readU(unsigned int & u) const {
    const bool retVal = (1 == fread(&u,sizeof(unsigned int),1,fp));
    if (!retVal) {
        ORSA_DEBUG("problems...");
    }
    return retVal;
}

bool RadioScienceGravityFile::readS(std::string & s) const {
    char line[8];
    const bool retVal = (8 == fread(&line,sizeof(char),8,fp));
    if (!retVal) {
        ORSA_DEBUG("problems...");
    }
    s = line;
    orsa::removeLeadingAndTrailingSpaces(s);
    return retVal;
}

void RadioScienceGravityFile::skipToNextRow() const {
    if ((ftell(fp)%RECORD_BYTES) != 0) {
        fseek(fp,RECORD_BYTES-(ftell(fp)%RECORD_BYTES),SEEK_CUR);
    }
}
