#include <orsaPDS/RadioScienceGravity.h>

#include <iostream>
#include <cstdio>

#include <orsa/print.h>

using namespace orsaPDS;

RadioScienceGravityFile::RadioScienceGravityFile(const std::string & fileName,
                                                 const size_t & RECORD_BYTES_,
                                                 const size_t & FILE_RECORDS_) :
    RECORD_BYTES(RECORD_BYTES_),
    FILE_RECORDS(FILE_RECORDS) {
    
    data = new RadioScienceGravityData;
    
    fp = fopen(fileName.c_str(),"r");
    if (!fp) {
        ORSA_DEBUG("cannot open file [%s]",fileName.c_str());
        return;
    }
    
    // read SHBDR_HEADER_TABLE
    {
        readD(data->R0);
        orsa::FromUnits(data->R0,orsa::Unit::KM);
        
        readD(data->GM);
        orsa::FromUnits(orsa::FromUnits(data->GM,orsa::Unit::KM,3),orsa::Unit::SECOND,-2);
        
        readD(data->sigmaGM);
        orsa::FromUnits(orsa::FromUnits(data->sigmaGM,orsa::Unit::KM,3),orsa::Unit::SECOND,-2);
        
        readI(data->degree);

        readI(data->order);

        readI(data->normalizationState);

        readI(data->numberOfCoefficients);
        
        readD(data->referenceLongitude);
        data->referenceLongitude *= orsa::degToRad();
        
        readD(data->referenceLatitude);
        data->referenceLatitude *= orsa::degToRad();
        
        skipToNextRow();
    }

    ORSA_DEBUG("n. coeff: %i",data->numberOfCoefficients);
    
    // read SHBDR_NAMES_TABLE
    {
        for (unsigned int k=0; k<data->numberOfCoefficients; ++k) {
            std::string s;
            readS(s);
            ORSA_DEBUG("name: [%s]",s.c_str());
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

bool RadioScienceGravityFile::readS(std::string & s) const {
    char line[8];
    const bool retVal = (8 == fread(&line,sizeof(char),8,fp));
    if (!retVal) {
        ORSA_DEBUG("problems...");
    }
    s = line;
    return retVal;
}

void RadioScienceGravityFile::skipToNextRow() const {
    long pos = ftell(fp);
    ORSA_DEBUG("pos before seek: %i",pos);
    if ((pos%RECORD_BYTES) != 0) {
        fseek(fp,RECORD_BYTES-(pos%RECORD_BYTES),SEEK_CUR);
    }
    pos = ftell(fp);
    ORSA_DEBUG("pos after seek: %i",ftell(fp));
}
