#include "simplex.h"

// #include "vesta.h"
// #include "gaskell.h"
// #include "eros_shape.h"

int main(int argc, char **argv) {
    
    const double km = orsa::FromUnits(1,orsa::Unit::KM);
    
    if (argc != 3) {
        printf("Usage: %s <vertex-file> <GeodesicGrid-Nsub>\n",argv[0]);
        exit(0);
    }
    
    const std::string inputFile = argv[1];
    const size_t Nsub = strtoul(argv[2],0,10);
    
    // geodesic grid
    orsa::TriShape::VertexVector v;
    orsa::TriShape::FaceVector f;
    orsa::TriShape::GeodesicGrid(v,f,Nsub);
    
    /* osg::ref_ptr<VestaShape> shapeModel = new VestaShape;
       if (!shapeModel->read(inputFile)) {
       ORSA_ERROR("problems encountered while reading shape file...");
       exit(0);
       }
    */
    
    /* osg::ref_ptr<ErosShape> shapeModel = new ErosShape;
       if (!shapeModel->read(inputFile)) {
       ORSA_ERROR("problems encountered while reading shape file...");
       exit(0);
       }
    */
    
    /*
    osg::ref_ptr<GaskellPlateModel> shapeModel = new GaskellPlateModel;
    if (!shapeModel->read(inputFile)) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    */
    
    orsa::TriShape::VertexVector inVec;
    {
        FILE * fp = fopen(inputFile.c_str(),"r");
        if (!fp) {
            ORSA_ERROR("problems encountered while reading vertex file [%s]",inputFile.c_str());
            exit(0);
        }
        double vx,vy,vz;
        char line[4096];
        while (fgets(line,4096,fp)) {
            if (3==sscanf(line,"%lf %lf %lf",&vx,&vy,&vz)) {
                inVec.push_back(km*orsa::Vector(vx,vy,vz));
            }
        }
        fclose(fp);
    }
    
    ORSA_DEBUG("GeodesicGrid v.size(): %u",v.size());
    ORSA_DEBUG("%s v.size(): %u",inputFile.c_str(),inVec.size());
    
    // process
    {
        // vars
        double      prod;
        double good_prod;
        orsa::Vector good_vec = orsa::Vector(0,0,1);
        for (unsigned int k=0; k<v.size(); ++k) {
            if (k%100==0) {
                ORSA_DEBUG("progress: %i/%i",k,v.size());
            }
            // init vars in loop
            good_prod = -1.0;
            good_vec  = orsa::Vector(0,0,1);
            for (unsigned int j=0; j<inVec.size(); ++j) {
                prod = v[k].normalized()*inVec[j].normalized();
                if (prod>good_prod) {
                    good_prod=prod;
                    good_vec =inVec[j];
                }
            }
            // save
            v[k] = good_vec;
        }
    }
    
    // write file in Gaskell format
    {
        char outFileName[4096];
        sprintf(outFileName,"%s_reShapeVec_%u.out",inputFile.c_str(),Nsub);
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
