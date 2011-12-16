#include "CubicChebyshevMassDistribution.h"

#include <orsa/unit.h>

std::vector< std::vector< std::vector<size_t> > > CubicChebyshevMassDistribution::indexTable;

size_t CubicChebyshevMassDistribution::totalSize(const size_t & degree) {
    return (degree+1)*(degree+2)*(degree+3)/6;
}

size_t CubicChebyshevMassDistribution::degree(const CoefficientType & coeff) {
    if (coeff.size() == 0) {
        ORSA_DEBUG("PROBLEM: degree not defined when coeff.size() == 0");
        return 0;
    } else {
        return coeff.size()-1;
    }
}

size_t CubicChebyshevMassDistribution::degree() const { 
    return degree(coeff);
}

void CubicChebyshevMassDistribution::resize(CoefficientType & coeff, const size_t & degree) {
    coeff.resize(degree+1);
    for (size_t i=0; i<=degree; ++i) {
        coeff[i].resize(degree+1-i);
        for (size_t j=0; j<=degree-i; ++j) {
            coeff[i][j].resize(degree+1-i-j);
        }
    }
}

size_t CubicChebyshevMassDistribution::index(const size_t & nx, const size_t & ny, const size_t & nz) {
    const size_t requestedDegree=nx+ny+nz;
    updateIndexTable(requestedDegree);
    return indexTable[nx][ny][nz];
}

void CubicChebyshevMassDistribution::triIndex(size_t & nx, size_t & ny, size_t & nz, const size_t & index) {
    size_t localDegree=0;
    while (CubicChebyshevMassDistribution::totalSize(localDegree)<=index) { ++localDegree; }
    updateIndexTable(localDegree);
    for (nx=0; nx<=localDegree; ++nx) {
        for (ny=0; ny<=localDegree; ++ny) {
            for (nz=0; nz<=localDegree; ++nz) {
                if (nx+ny+nz==localDegree) {
                    // ORSA_DEBUG("nx: %i  ny: %i  nz: %i  index: %i  localDegree: %i",nx,ny,nz,index,localDegree);
                    if (indexTable[nx][ny][nz]==index) return;
                }
            }
        }
    }
    ORSA_DEBUG("index [%i] not found",index);
}

void CubicChebyshevMassDistribution::updateIndexTable(const size_t & requestedDegree) {
    if (indexTable.size() < (requestedDegree+1)) {
        indexTable.resize(requestedDegree+1);
        for (size_t i=0; i<=requestedDegree; ++i) {
            indexTable[i].resize(requestedDegree+1-i);
            for (size_t j=0; j<=requestedDegree-i; ++j) {
                indexTable[i][j].resize(requestedDegree+1-i-j);
            }
        }
        size_t idx=0;
        size_t localDegree=0;
        while (localDegree <= requestedDegree) {
            for (unsigned int i=0; i<=localDegree; ++i) {
                for (unsigned int j=0; j<=localDegree; ++j) {
                    for (unsigned int k=0; k<=localDegree; ++k) {
                        if (i+j+k==localDegree) {
                            // ORSA_DEBUG("inserting %i-%i-%i  index: %i",i,j,k,idx);
                            indexTable[i][j][k] = idx++;
                        }
                    }
                }
            }
            ++localDegree;
        }
    }
}

CubicChebyshevMassDistribution::CubicChebyshevMassDistribution(const CoefficientType & coefficient,
                                                               const double & densityScale_,
                                                               const double & R0,
                                                               const LayerData * layerData_) :
    orsa::MassDistribution(),
    coeff(coefficient),
    densityScale(densityScale_),
    oneOverR0(1.0/R0),
    layerData(layerData_) {
}

CubicChebyshevMassDistribution::~CubicChebyshevMassDistribution() { }

double CubicChebyshevMassDistribution::density(const orsa::Vector & p) const {
    if (0) {
        // debug
        static size_t calls=0;
        ++calls;
        if (calls%1000==0) ORSA_DEBUG("calls: %i",calls);
    }
    // if (coeff.size() == 0) return 0.0;
    if ( (coeff.size() == 0) &&
         (layerData.get() == 0) ) {
        return 0.0;
    }
    std::vector<double> Tx, Ty, Tz;
    const double Tx_arg = p.getX()*oneOverR0;
    const double Ty_arg = p.getY()*oneOverR0;
    const double Tz_arg = p.getZ()*oneOverR0;
    //
    if ( (fabs(Tx_arg) > 1.0) ||
         (fabs(Ty_arg) > 1.0) ||
         (fabs(Tz_arg) > 1.0) ) {
        ORSA_DEBUG("problem: R0 value too small");
        exit(0);
    }
    //
    const size_t T_degree = degree();
    orsa::ChebyshevT(Tx,T_degree,Tx_arg);
    orsa::ChebyshevT(Ty,T_degree,Ty_arg);
    orsa::ChebyshevT(Tz,T_degree,Tz_arg);
    double density = 0.0;
    for (size_t i=0; i<=T_degree; ++i) {
        for (size_t j=0; j<=T_degree-i; ++j) {
            for (size_t k=0; k<=T_degree-i-j; ++k) {
                density += coeff[i][j][k]*Tx[i]*Ty[j]*Tz[k];
            }
        }
    }
    density *= densityScale;
    
    if (layerData.get() != 0) {
        density += layerData->density(p);
    }
    
    return density;
}

/***/


// decompose a generic mass distribution into a cubic chebyshev mass distribuiton
// note: there are no tests on whether the points tested are inside or outside the body shape,
//       but that should not matter; what matters is that the new mass distribution
//       returns the same density as the input at any given point
CubicChebyshevMassDistribution * CubicChebyshevMassDistributionDecomposition(const orsa::MassDistribution * massDistribution,
                                                                             const size_t & degree,
                                                                             const double & densityScale,
                                                                             const double & R0,
                                                                             const LayerData * layerData) {
    CubicChebyshevMassDistribution::CoefficientType coeff;
    CubicChebyshevMassDistribution::resize(coeff,degree);
    
    for (size_t running_n=0; running_n<=degree; ++running_n) {
        ORSA_DEBUG("[%zi/%zi]",running_n,degree);
        for (size_t nx=0; nx<=degree; ++nx) {
            for (size_t ny=0; ny<=degree-nx; ++ny) {
                for (size_t nz=0; nz<=degree-nx-ny; ++nz) {
                    if (nx+ny+nz != running_n) continue;
                    double sum = 0.0;
                    for (size_t kx=0; kx<degree; ++kx) {
                        const double x = cos(orsa::pi()*(kx+0.5)/degree);
                        for (size_t ky=0; ky<degree; ++ky) {
                            const double y = cos(orsa::pi()*(ky+0.5)/degree);
                            for (size_t kz=0; kz<degree; ++kz) {
                                const double z = cos(orsa::pi()*(kz+0.5)/degree);
                                
                                /* sum +=
                                   cos(nx*orsa::pi()*(kx+0.5)/degree) *
                                   cos(ny*orsa::pi()*(ky+0.5)/degree) *
                                   cos(nz*orsa::pi()*(kz+0.5)/degree) *
                                   massDistribution->density(R0*orsa::Vector(x,y,z)) / densityScale;
                                */
                                
                                double delta_density = 
                                    massDistribution->density(R0*orsa::Vector(x,y,z));
                                if (layerData != 0) {
                                    delta_density -=
                                        layerData->density(R0*orsa::Vector(x,y,z));
                                }
                                
                                sum +=
                                    cos(nx*orsa::pi()*(kx+0.5)/degree) *
                                    cos(ny*orsa::pi()*(ky+0.5)/degree) *
                                    cos(nz*orsa::pi()*(kz+0.5)/degree) *
                                    delta_density / densityScale;
                            }
                        }
                    }
                    coeff[nx][ny][nz] =
                        (2-orsa::kronecker(0,nx)) *
                        (2-orsa::kronecker(0,ny)) *
                        (2-orsa::kronecker(0,nz)) *
                        sum/orsa::cube(degree);
                    /* ORSA_DEBUG("coeff[%02i][%02i][%02i] = %20.12f",
                       nx,ny,nz,coeff[nx][ny][nz]);
                    */
                }
            }
        }
    }
    
    CubicChebyshevMassDistribution * CCMD =
        new CubicChebyshevMassDistribution(coeff,densityScale,R0,layerData);
    
    return CCMD;
}

/***/

bool CubicChebyshevMassDistributionFile::read(CubicChebyshevMassDistributionFile::DataContainer & data, const std::string & fileName) {
    FILE * fp = fopen(fileName.c_str(),"r");
    if (!fp) {
        ORSA_DEBUG("cannot open file [%s]",fileName.c_str());
        return false;
    }
    data.clear();
    CubicChebyshevMassDistributionFile::DataType dataElement;
    while (read(dataElement,fp)) {
        data.push_back(dataElement);
    }
    fclose(fp);
    return true;
}

bool CubicChebyshevMassDistributionFile::read(CubicChebyshevMassDistributionFile::DataContainer & data, const std::string & fileName, const double & limitDeltaDensity) {
    FILE * fp = fopen(fileName.c_str(),"r");
    if (!fp) {
        ORSA_DEBUG("cannot open file [%s]",fileName.c_str());
        return false;
    }
    data.clear();
    size_t num_read=0, num_imported=0;
    CubicChebyshevMassDistributionFile::DataType dataElement;
    while (read(dataElement,fp)) {
        ++num_read;
        if (dataElement.deltaDensity <= limitDeltaDensity) {
            ++num_imported;
            data.push_back(dataElement);
        }
    }
    fclose(fp);
    if (num_read != 0) ORSA_DEBUG("imported: %i/%i [%.3f%%]",num_imported,num_read,100*(double)num_imported/(double)num_read);
    return true;
}

bool CubicChebyshevMassDistributionFile::write(const CubicChebyshevMassDistributionFile::DataContainer & data, const std::string & fileName) {
    FILE * fp = fopen(fileName.c_str(),"w");
    if (!fp) {
        ORSA_DEBUG("cannot open file [%s]",fileName.c_str());
        return false;
    }
    CubicChebyshevMassDistributionFile::DataContainer::const_iterator it = data.begin();
    while (it != data.end()) {
        write((*it),fp);
        ++it;
    }
    fclose(fp);
    return true;
}

bool CubicChebyshevMassDistributionFile::write(const CubicChebyshevMassDistributionFile::DataType & data, const std::string & fileName) {
    FILE * fp = fopen(fileName.c_str(),"w");
    if (!fp) {
        ORSA_DEBUG("cannot open file [%s]",fileName.c_str());
        return false;
    }
    write(data,fp);
    fclose(fp);
    return true;
}

bool CubicChebyshevMassDistributionFile::append(const CubicChebyshevMassDistributionFile::DataType & data, const std::string & fileName) {
    FILE * fp = fopen(fileName.c_str(),"a");
    if (!fp) {
        ORSA_DEBUG("cannot open file [%s]",fileName.c_str());
        return false;
    }
    write(data,fp);
    fclose(fp);
    return true;
}

bool CubicChebyshevMassDistributionFile::read(CubicChebyshevMassDistributionFile::DataType & data, FILE * fp) {
    if (1 != gmp_fscanf(fp,"%lf",&data.minDensity)) return false;                                
    data.minDensity = orsa::FromUnits(orsa::FromUnits(data.minDensity,orsa::Unit::GRAM),orsa::Unit::CM,-3);
    if (1 != gmp_fscanf(fp,"%lf",&data.maxDensity)) return false;                                
    data.maxDensity = orsa::FromUnits(orsa::FromUnits(data.maxDensity,orsa::Unit::GRAM),orsa::Unit::CM,-3);
    if (1 != gmp_fscanf(fp,"%lf",&data.deltaDensity)) return false;                                    
    data.deltaDensity = orsa::FromUnits(orsa::FromUnits(data.deltaDensity,orsa::Unit::GRAM),orsa::Unit::CM,-3);
    if (1 != gmp_fscanf(fp,"%lf",&data.penalty)) return false;
    if (1 != gmp_fscanf(fp,"%lf",&data.densityScale)) return false;                                
    data.densityScale = orsa::FromUnits(orsa::FromUnits(data.densityScale,orsa::Unit::GRAM),orsa::Unit::CM,-3);
    if (1 != gmp_fscanf(fp,"%lf",&data.R0)) return false;                                
    data.R0 = orsa::FromUnits(data.R0,orsa::Unit::KM,1);
    if (1 != gmp_fscanf(fp,"%zi",&data.SH_degree)) return false;   
    size_t T_degree;
    if (1 != gmp_fscanf(fp,"%zi",&T_degree)) return false;                                
    CubicChebyshevMassDistribution::resize(data.coeff,T_degree);
    for (size_t runningDegree=0; runningDegree<=T_degree; ++runningDegree) {
        for (size_t i=0; i<=T_degree; ++i) {
            for (size_t j=0; j<=T_degree-i; ++j) {
                for (size_t k=0; k<=T_degree-i-j; ++k) {
                    if (i+j+k == runningDegree) {
                        if (1 != gmp_fscanf(fp,"%lf",&data.coeff[i][j][k])) return false;
                    }
                }
            }
        }
    }
    bool haveLayerData=false;
    /* double baseDensity;
       if (1 == gmp_fscanf(fp,"%lf",&baseDensity)) {
       haveLayerData=true;
       baseDensity = orsa::FromUnits(orsa::FromUnits(baseDensity,orsa::Unit::GRAM),orsa::Unit::CM,-3);
       }
    */
    LayerData::EllipsoidLayerVectorType ellipsoidLayerVector;
    size_t ellipsoidLayerVectorSize;
    if (1 == gmp_fscanf(fp,"%zi",&ellipsoidLayerVectorSize)) {
        haveLayerData=true;
        double excessDensity;
        double a,b,c;
        double v0x,v0y,v0z;
        for (unsigned int k=0; k<ellipsoidLayerVectorSize; ++k) {
            if (7 == gmp_fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",
                                &excessDensity,
                                &a,
                                &b,
                                &c,
                                &v0x,
                                &v0y,
                                &v0z)) {
                excessDensity = orsa::FromUnits(orsa::FromUnits(excessDensity,orsa::Unit::GRAM),orsa::Unit::CM,-3);
                a = orsa::FromUnits(a,orsa::Unit::KM);
                b = orsa::FromUnits(b,orsa::Unit::KM);
                c = orsa::FromUnits(c,orsa::Unit::KM);
                v0x = orsa::FromUnits(v0x,orsa::Unit::KM);
                v0y = orsa::FromUnits(v0y,orsa::Unit::KM);
                v0z = orsa::FromUnits(v0z,orsa::Unit::KM);
                ellipsoidLayerVector.push_back(new LayerData::EllipsoidLayer(excessDensity,a,b,c,orsa::Vector(v0x,v0y,v0z)));
            } else {
                return false;
            }
        }
    }
    if (haveLayerData) {
        // data.layerData = new LayerData(baseDensity,ellipsoidLayerVector);
        data.layerData = new LayerData(ellipsoidLayerVector);
    }
    return true;
}

bool CubicChebyshevMassDistributionFile::write(const CubicChebyshevMassDistributionFile::DataType & data, FILE * fp) {
    gmp_fprintf(fp,"%.3f ",orsa::FromUnits(orsa::FromUnits(data.minDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
    gmp_fprintf(fp,"%.3f ",orsa::FromUnits(orsa::FromUnits(data.maxDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
    gmp_fprintf(fp,"%.3f ",orsa::FromUnits(orsa::FromUnits(data.deltaDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
    gmp_fprintf(fp,"%.3f ",data.penalty);
    gmp_fprintf(fp,"%.3f ",orsa::FromUnits(orsa::FromUnits(data.densityScale,orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
    gmp_fprintf(fp,"%g ",orsa::FromUnits(data.R0,orsa::Unit::KM,-1));
    gmp_fprintf(fp,"%i ",data.SH_degree);
    const size_t T_degree = CubicChebyshevMassDistribution::degree(data.coeff);
    gmp_fprintf(fp,"%i ",T_degree);
    for (size_t runningDegree=0; runningDegree<=T_degree; ++runningDegree) {
        for (size_t i=0; i<=T_degree; ++i) {
            for (size_t j=0; j<=T_degree-i; ++j) {
                for (size_t k=0; k<=T_degree-i-j; ++k) {
                    if (i+j+k == runningDegree) {
                        gmp_fprintf(fp,"%+9.6f ",data.coeff[i][j][k]);
                    }
                }
            }
        }
    }
    if (data.layerData.get() != 0) {
        // gmp_fprintf(fp,"%.3f ",orsa::FromUnits(orsa::FromUnits(data.layerData->baseDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
        const LayerData::EllipsoidLayerVectorType & lv = data.layerData->ellipsoidLayerVector;
        if (lv.size() > 0) {
            gmp_fprintf(fp,"%i ",lv.size());
            for (unsigned int k=0; k<lv.size(); ++k) {
                gmp_fprintf(fp,"%.3f %.g %.g %.g %.g %.g %.g ",
                            orsa::FromUnits(orsa::FromUnits(lv[k]->excessDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
                            orsa::FromUnits(lv[k]->a,orsa::Unit::KM,-1),
                            orsa::FromUnits(lv[k]->b,orsa::Unit::KM,-1),
                            orsa::FromUnits(lv[k]->c,orsa::Unit::KM,-1),
                            orsa::FromUnits(lv[k]->v0.getX(),orsa::Unit::KM,-1),
                            orsa::FromUnits(lv[k]->v0.getY(),orsa::Unit::KM,-1),
                            orsa::FromUnits(lv[k]->v0.getZ(),orsa::Unit::KM,-1));
            }
        }
    }
    gmp_fprintf(fp,"\n");
    return true;
}
