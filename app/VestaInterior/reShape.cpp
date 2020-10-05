#include "simplex.h"

#include "shape.h"

int main(int argc, char **argv) {
    
    if (argc != 3) {
        printf("Usage: %s <plate-model-file-in> <GeodesicGrid-Nsub>\n",argv[0]);
        exit(0);
    }
    
    const std::string inputFile = argv[1];
    const size_t Nsub = strtoul(argv[2],0,10);
    
    // geodesic grid
    orsa::TriShape::VertexVector v;
    orsa::TriShape::FaceVector f;
    orsa::TriShape::GeodesicGrid(v,f,Nsub);
    
    // PICK one input file type
    
    
    osg::ref_ptr<InputShape> shapeModel = new InputShape;
       if (!shapeModel->read(inputFile)) {
       ORSA_ERROR("problems encountered while reading shape file...");
       exit(0);
    }
        
    // *****
    
    ORSA_DEBUG("GeodesicGrid v.size(): %u",v.size());
    ORSA_DEBUG("%s v.size(): %u",inputFile.c_str(),shapeModel->getVertexVector().size());
    
    // process
    {
        orsa::Vector intersectionPoint;
        orsa::Vector normal;
        const orsa::Vector P(0,0,0);
        const bool fullLine = false;
        for (unsigned int k=0; k<v.size(); ++k) {
            if (k%100==0) {
                ORSA_DEBUG("progress: %i/%i",k,v.size());
            }
            if (!shapeModel->rayIntersection(intersectionPoint,normal,P,v[k],fullLine)) {
                ORSA_DEBUG("problems...");
                exit(0);
            } else {
                v[k] = intersectionPoint - P;
            }
        }
    }
    
    // write file in Gaskell format
    {
        char outFileName[4096];
        sprintf(outFileName,"%s_reShape_%u.out",inputFile.c_str(),Nsub);
        FILE * fp = fopen(outFileName,"w");
        ORSA_DEBUG("writing file [%s]...",outFileName);
        fprintf(fp,"%i\n",v.size());
        for (size_t j=0; j<v.size(); ++j) {
            fprintf(fp,"%i %10.5f %10.5f %10.5f\n",
                    j+1,
                    orsa::FromUnits(v[j].getX(),orsa::Unit::KM,-1),
                    orsa::FromUnits(v[j].getY(),orsa::Unit::KM,-1),
                    orsa::FromUnits(v[j].getZ(),orsa::Unit::KM,-1));
        }
        fprintf(fp,"%i\n",f.size());
        for (size_t p=0; p<f.size(); ++p) {
            fprintf(fp,"%i %i %i %i\n",
                    1+p,
                    1+f[p].i(),
                    1+f[p].j(),
                    1+f[p].k());
        }
        fclose(fp);
    }
    
    return 0;
}
