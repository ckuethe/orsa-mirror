#include "simplex.h"

// #include "vesta.h"
// #include "gaskell.h"
// #include "eros_shape.h"
#include "shape.h"

#include "mpreal.h"

template <typename T> std::vector< std::vector< std::vector<size_t> > > SimplexIntegration<T>::indexTable;
template <typename T> std::vector< std::vector< std::vector< std::vector<size_t> > > > SimplexIntegration<T>::index4Table;

int main(int argc, char **argv) {
    
    typedef mpfr::mpreal F;
    mpfr::mpreal::set_default_prec(128); // number of bits of precision
    
    if (argc != 7) {
        printf("Usage: %s <plate-model-file> <R0_km> <min-degree> <max-degree> <mod-N> <mod-i>\n",argv[0]);
        exit(0);
    }
    
    const std::string inputFile = argv[1];
    const double R0 = orsa::FromUnits(atof(argv[2]),orsa::Unit::KM);
    const int minDegree = atoi(argv[3]);
    const int maxDegree = atoi(argv[4]);
    const int mod_N = atoi(argv[5]);
    const int mod_i = atoi(argv[6]);
    
    if ( (R0 <= 0.0) ||
         (minDegree < 0) ||
         (maxDegree < minDegree) ||
         (mod_N < 1) ||
         (mod_i < 0) ||
         (mod_i >= mod_N) ) {
        ORSA_DEBUG("invalid input...");
        exit(0);
    }

    // safer over NFS
    // sqlite3_vfs_register(sqlite3_vfs_find("unix-dotfile"), 1);
    
    const std::string SQLiteDBFileName = getSqliteDBFileName_simplex(inputFile,R0);
	ORSA_DEBUG("SQLiteDBFileName: [%s]",SQLiteDBFileName.c_str());
    
    orsa::Debug::instance()->initTimer();
         
    osg::ref_ptr<InputShape> shapeModel = new InputShape;
       if (!shapeModel->read(inputFile)) {
       ORSA_ERROR("problems encountered while reading shape file...");
       exit(0);
    }
    
    if (1) {
        // output .xyz file for plotting with GMT
        char filename[1024];
        sprintf(filename,"%s.xyz",inputFile.c_str());
        FILE * fp;
        fp = fopen(filename,"r");
        if (fp != 0) {
            ORSA_DEBUG("file [%s] already existing, skipping",filename);
            fclose(fp);
        } else {
            fp = fopen(filename,"w");
            ORSA_DEBUG("writing file [%s]",filename);
            orsa::TriShape::VertexVector vertex = shapeModel->getVertexVector();
            orsa::TriShape::VertexVector::const_iterator it = vertex.begin();
            while (it != vertex.end()) {
                const orsa::Vector & v = (*it);
                const double lat = asin(v.getZ()/v.length());
                const double lon = fmod(orsa::twopi()+atan2(v.getY(),v.getX()),orsa::twopi());
                fprintf(fp,"%g %g %g\n",
                        lon*orsa::radToDeg(),
                        lat*orsa::radToDeg(),
                        orsa::FromUnits(v.length(),orsa::Unit::KM,-1));
                ++it;
            }
            fclose(fp);
        }
    }
    
    osg::ref_ptr<SimplexIntegration<F> > si = new SimplexIntegration<F>(shapeModel.get(), R0, SQLiteDBFileName);
    // osg::ref_ptr<SimplexIntegration> si_unit_R0 = new SimplexIntegration(vestaShape.get(),1.0);
    
    si->reserve(maxDegree);
    // si_unit_R0->reserve(maxDegree);
    for (size_t degree=minDegree; degree<=(size_t)maxDegree; ++degree) {
        for (size_t i=0; i<=degree; ++i) {
            for (size_t j=0; j<=degree; ++j) {
                for (size_t k=0; k<=degree; ++k) {
                    if (i+j+k==degree) {
                        const size_t index = SimplexIntegration<F>::getIndex(i,j,k);
                        if ((index%mod_N)==(size_t)mod_i) {

                            const double integral_ijk = si->getIntegral(i,j,k);
                            ORSA_DEBUG("integral [%02i,%02i,%02i]: %+20.12f",i,j,k,integral_ijk);
                            // ORSA_DEBUG("scaled: %+16.6e",integral_ijk*orsa::int_pow(R0,3+degree)); // 3=jacobian, degree=transformation
                        }
                    }
                }
            }
        }
    }
    
    return 0;
}
