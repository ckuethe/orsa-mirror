#include "simplex.h"
#include "SH2ijk.h"

#include <orsa/legendre.h>
#include <orsa/util.h>

// norm_coeff = normalization_factor * coeff
double normalization_factor(const size_t & l,
                            const size_t & m) {
    return orsa::normalization_sphericalHarmonicsToNormalizedSphericalHarmonics(l,m).get_d();
}

double radius(const double & theta,
              const double & phi,
              const SHcoeff & norm_A,
              const SHcoeff & norm_B) {
    const double c_theta = cos(theta);
    double r=0;
    if (norm_A.size() != norm_B.size()) {
        ORSA_DEBUG("problems...");
    }
    for (size_t l=0; l<norm_A.size(); ++l) {
        if (norm_A[l].size() != norm_B[l].size()) {
            ORSA_DEBUG("problems...");
        }
        for (size_t m=0; m<norm_A[l].size(); ++m) {
            if ( (norm_A[l][m] != 0.0) ||
                 (norm_B[l][m] != 0.0) ) {
                const double Nf = normalization_factor(l,m);
                const double Plm = orsa::LegendreP(l,m,c_theta).get_d()/Nf;
                // ORSA_DEBUG("LegendreP[%i][%i](%g) = %g",l,m,c_theta,Plm);
                r += Plm*(norm_A[l][m]*cos(m*phi));
                if (m!=0) {
                    r += Plm*(norm_B[l][m]*sin(m*phi));
                }   
            }
        }
    }
    return r;
}

double radius(orsa::Vector & r, const SHcoeff & norm_A, const SHcoeff & norm_B) {
    const orsa::Vector u = r.normalized();
    const double theta = acos(u.getZ());
    const double phi   = atan2(u.getY(),u.getX());
    //
    r = u*radius(theta,phi,norm_A,norm_B);
    return r.length();
}

int main(int argc, char **argv) {
    
    const double km = orsa::FromUnits(1,orsa::Unit::KM);
    
    if (argc != 4) {
        printf("Usage: %s <SH-file> <max-degree> <GeodesicGrid-Nsub>\n",argv[0]);
        exit(0);
    }
    
    const std::string inputFile = argv[1];
    const size_t max_degree = strtoul(argv[2],0,10);
    const size_t Nsub       = strtoul(argv[3],0,10);
    
    // geodesic grid
    orsa::TriShape::VertexVector v;
    orsa::TriShape::FaceVector f;
    orsa::TriShape::GeodesicGrid(v,f,Nsub);
    
    SHcoeff norm_A;
    SHcoeff norm_B;
    readSH(norm_A,norm_B,inputFile,orsa::Unit::KM);
    if (norm_A.size()>(max_degree+1)) norm_A.resize(max_degree+1);
    if (norm_B.size()>(max_degree+1)) norm_B.resize(max_degree+1);
    
    /*
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
    */
    
    ORSA_DEBUG("GeodesicGrid v.size(): %u",v.size());
    // ORSA_DEBUG("%s v.size(): %u",inputFile.c_str(),inVec.size());
    
    // set vectors
    for (unsigned int k=0; k<v.size(); ++k) {
        radius(v[k],norm_A,norm_B);
    }
    
    // write file in Gaskell format
    {
        
        char outFileName[4096];
        sprintf(outFileName,"%s_trim_%u_reShapeSH_%u.out",inputFile.c_str(),max_degree,Nsub);
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
