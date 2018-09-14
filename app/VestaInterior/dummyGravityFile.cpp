#include <orsaPDS/RadioScienceGravity.h>
#include <orsa/unit.h>

int main(int argc, char **argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if (argc != 3) {
        printf("Usage: %s <R0_km> <degree>\n",argv[0]);
        exit(0);
    }
    
    const double R0     = orsa::FromUnits(atof(argv[1]),orsa::Unit::KM);
    const size_t degree = atoi(argv[2]);
    
    osg::ref_ptr<orsaPDS::RadioScienceGravityData> gravityData = new orsaPDS::RadioScienceGravityData;
    gravityData->R0     = R0;
    gravityData->degree = degree;
    gravityData->order  = degree;
    gravityData->numberOfCoefficients = (degree+1)*(degree+1)-3;
    
    {
        unsigned int index = 0;
        gravityData->hash["GM"] = index++;
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
    
    orsaPDS::RadioScienceGravityFile::write(gravityData.get(),"dummyGravityFile.bin",512,1518);
    
    return 0;
}
