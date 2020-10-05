#include <orsaPDS/RadioScienceGravity.h>
#include <orsa/unit.h>

int main(int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if (argc != 5) {
        printf("Usage: %s <input-file> <R0_km> <GM [km^3/s^2]> <degree>\n",argv[0]);
        exit(0);
    }
    
    const std::string inputFile = argv[1];
    const double R0     = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const double GM     = orsa::FromUnits(orsa::FromUnits(atof(argv[3]),orsa::Unit::KM,3),orsa::Unit::SECOND,-2);
    const int    degree = atoi(argv[4]);
    
    osg::ref_ptr<orsaPDS::RadioScienceGravityData> gravityData = new orsaPDS::RadioScienceGravityData;
    gravityData->R0     = R0;
    gravityData->GM     = GM;
    gravityData->degree = degree;
    gravityData->order  = degree;
    gravityData->numberOfCoefficients = (degree+1)*(degree+1)-3;
    
    {
        unsigned int index = 0;
        gravityData->hash["GM"] = index++;
        // gravityData->hash[orsaPDS::RadioScienceGravityData::keyC(0,0)] = index++; // C00 for uncertainty of GM
        for (unsigned int l=2; l<=gravityData->degree; ++l) {
            for (unsigned int m=0; m<=l; ++m) {
                gravityData->hash[orsaPDS::RadioScienceGravityData::keyC(l,m)] = index++;
                if (m != 0) {
                    gravityData->hash[orsaPDS::RadioScienceGravityData::keyS(l,m)] = index++;
                }       
            }
        }        
    }
    
    gravityData->coeff.resize(gravityData->numberOfCoefficients);
    for (unsigned int j=0; j<gravityData->numberOfCoefficients; ++j) {
        gravityData->coeff[j] = 0.0;
    }
    
    gravityData->covar.resize(gravityData->numberOfCoefficients);
    for (unsigned int j=0; j<gravityData->numberOfCoefficients; ++j) {
        gravityData->covar[j].resize(j+1);
        for (unsigned int k=0; k<=j; ++k) {
            gravityData->covar[j][k] = 0.0;
        }
    }

    { 
        // now finally read gravity data from file

        FILE * fp = fopen(inputFile.c_str(),"r");
        if (!fp) {
            ORSA_DEBUG("problems: cannot open file [%s]",inputFile.c_str());
            exit(0);
        }

        char line[4096];
        int l, m;
        double Clm, Slm, sigmaClm, sigmaSlm;
        while (fgets(line,4096,fp)) {
            // ORSA_DEBUG("line: [%s]",line);
            if (6 == sscanf(line,"%i %i %lf %lf %lf %lf",&l,&m,&Clm,&Slm,&sigmaClm,&sigmaSlm)) {
                // normalizations = ?
                if ((l==0) && (m==0)) {
                    // for uncertainty of GM
                    // gravityData->setCovar(orsaPDS::RadioScienceGravityData::keyC(l,m),orsaPDS::RadioScienceGravityData::keyC(l,m),orsa::square(sigmaClm));
                    // gravityData->setCovar("GM","GM",orsa::square(sigmaClm));
                    gravityData->setCovar("GM","GM",orsa::square(sigmaClm*GM/Clm)); // include Clm and GM for scaling...
                }
                if ((l>=2) &&
                    (l<=degree)) {
                    gravityData->setCoeff(orsaPDS::RadioScienceGravityData::keyC(l,m),Clm);
                    gravityData->setCovar(orsaPDS::RadioScienceGravityData::keyC(l,m),orsaPDS::RadioScienceGravityData::keyC(l,m),orsa::square(sigmaClm));
                    if (m != 0) {  
                        gravityData->setCoeff(orsaPDS::RadioScienceGravityData::keyS(l,m),Slm);
                        gravityData->setCovar(orsaPDS::RadioScienceGravityData::keyS(l,m),orsaPDS::RadioScienceGravityData::keyS(l,m),orsa::square(sigmaSlm));
                    }
                }
            }
        }
        
        fclose(fp);
    }
    
    orsaPDS::RadioScienceGravityFile::write(gravityData.get(),"gravityFromAsciiFile.bin",512,1518);
    
    return 0;
}
