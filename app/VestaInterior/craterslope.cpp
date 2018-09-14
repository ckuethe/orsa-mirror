#include <orsa/debug.h>
#include <orsa/print.h>
#include <orsa/unit.h>
#include <orsa/statistic.h>

#include "gaskell.h"

class vertex {
public:
    orsa::Cache<orsa::Vector> v, u;
    orsa::Cache<double> lon, lat;
    orsa::Cache<double> slope;
public:
    void print() const {
        orsa::print(v);
        orsa::print(lon);
        orsa::print(lat);
        orsa::print(slope);
    }
};

class crater {
public:
    orsa::Cache<double> lon, lat, diam;
public:
    void print() const {
        orsa::print(lon);
        orsa::print(lat);
        orsa::print(diam);
    }
};

int main (int argc, char **argv) {
    
    if (argc != 3) {
        ORSA_DEBUG("Usage: %s <slope_support.dat> <crater-db-file>",argv[0]);
        exit(0);
    }
    
    const std::string slopeSupportFile = argv[1];
    const std::string craterFile = argv[2];
    
    std::vector<vertex> vertexVector;
    
    char line[1024];
    
    FILE * fp_slopeSupport = fopen(slopeSupportFile.c_str(),"r");
    if (fp_slopeSupport == 0) {
        ORSA_DEBUG("cannot open file [%s]",slopeSupportFile.c_str());
        exit(0);
    }
    
    {
        double Rx, Ry, Rz, lon, lat, slope;
        vertex vtx;
        while (fgets(line,1024,fp_slopeSupport) != 0) {
            if (6 == sscanf(line,"%lf %lf %lf %*s   %*s %*s %*s %*s   %lf %lf   %*s %*s %*s   %*s %*s %*s %*s %*s %*s %*s %*s   %*s %*s   %*s %*s %*s   %*s %lf %*s   %*s %*s %*s   %*s %*s %*s",
                            &Rx,&Ry,&Rz,&lon,&lat,&slope)) {
                Rx = orsa::FromUnits(Rx,orsa::Unit::KM);
                Ry = orsa::FromUnits(Ry,orsa::Unit::KM);
                Rz = orsa::FromUnits(Rz,orsa::Unit::KM);
                lon *= orsa::degToRad();
                lat *= orsa::degToRad();
                slope = (180.0-slope)*orsa::degToRad();
                
                vtx.v = orsa::Vector(Rx,Ry,Rz);
                vtx.u = (*vtx.v).normalized();
                vtx.lon = lon;
                vtx.lat = lat;
                vtx.slope = slope;
                vertexVector.push_back(vtx);
            }   
        }
    }
    
    fclose(fp_slopeSupport);
    
    std::vector<crater> craterVector;
    
    FILE * fp_crater = fopen(craterFile.c_str(),"r");
    if (fp_crater == 0) {
        ORSA_DEBUG("cannot open file [%s]",craterFile.c_str());
        exit(0);
    }

    {
        double lon, lat, diam;
        crater ctr;
        while (fgets(line,1024,fp_crater) != 0) {
            if (strlen(line) < 130) continue;
            if (line[0] == '#') continue;
            if (3 == sscanf(line,"%*s %*s %*s %*s %*s %*s %*s %*s %*s %lf %lf %lf",&lon,&lat,&diam)) {
                lon *= orsa::degToRad();
                lat *= orsa::degToRad();
                diam = orsa::FromUnits(diam,orsa::Unit::KM);
                
                ctr.lon = lon;
                ctr.lat = lat;
                ctr.diam = diam;
                craterVector.push_back(ctr);
            }
        }
    }
    
    fclose(fp_crater);

    /* ORSA_DEBUG("vertexVector.size(): %i   craterVector.size(): %i",
       vertexVector.size(),craterVector.size());
    */
    
    {
        double s_lon, c_lon;
        double s_lat, c_lat;
        osg::ref_ptr< orsa::Statistic<double> > stat_slope = new orsa::Statistic<double>;
        osg::ref_ptr< orsa::Statistic<double> > stat_R = new orsa::Statistic<double>;
        for (size_t ck=0; ck<craterVector.size(); ++ck) {
            const crater & ctr = craterVector[ck];
            // sincos(ctr.lon,&s_lon,&c_lon);
            s_lon = sin(ctr.lon);
			c_lon = cos(ctr.lon);
			// sincos(ctr.lat,&s_lat,&c_lat);
            s_lat = sin(ctr.lat);
			c_lat = cos(ctr.lat);
			const orsa::Vector uc(c_lon*c_lat,s_lon*c_lat,s_lat);
            vertex closestVertex = vertexVector[0];
            for (size_t vk=0; vk<vertexVector.size(); ++vk) {
                if ((vertexVector[vk].u*uc) > (closestVertex.u*uc)) {
                    closestVertex = vertexVector[vk];
                }
            }
            const double maxAngle = atan2(ctr.diam,2*(*closestVertex.v).length());
            const double minScalarProduct = cos(maxAngle);
            
            stat_slope->reset();
            stat_R->reset();
            for (size_t vk=0; vk<vertexVector.size(); ++vk) {
                if ((vertexVector[vk].u*uc) > minScalarProduct) {
                    stat_slope->insert(vertexVector[vk].slope);
                    stat_R->insert((*vertexVector[vk].v).length());
                }
            }
            
            // at least two points...
            if (stat_slope->entries() > 1) {
                ORSA_DEBUG("%g %g %g %g %g %g %g %g %g %g %g %g %Zi",
                           orsa::radToDeg()*ctr.lon,
                           orsa::radToDeg()*ctr.lat,
                           orsa::FromUnits(ctr.diam,orsa::Unit::KM,-1),
                           orsa::radToDeg()*stat_slope->min(),
                           orsa::radToDeg()*stat_slope->max(),
                           orsa::radToDeg()*(stat_slope->max()-stat_slope->min()),
                           orsa::radToDeg()*stat_slope->average(),
                           orsa::radToDeg()*(stat_slope->average()-stat_slope->min()),
                           orsa::radToDeg()*(stat_slope->max()-stat_slope->average()),
                           orsa::FromUnits(stat_R->average(),orsa::Unit::KM,-1),
                           orsa::FromUnits(stat_R->max()-stat_R->min(),orsa::Unit::KM,-1),
                           (stat_R->max()-stat_R->min())/ctr.diam,
                           stat_slope->entries().get_mpz_t());
            }
        }
    }
    
    return 0;
}
