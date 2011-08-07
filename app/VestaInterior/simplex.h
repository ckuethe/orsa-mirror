#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/shape.h>

#warning if shape is strongly concave and a simplex covers volume outside the body shape, then the results are incorrect (including volume computations...)

#warning default origin for 4th simplex vertex, should be a parameter of the class??                

class SimplexIntegration : public osg::Referenced {
public:
    SimplexIntegration(const orsa::TriShape * s,
                       const double & R0_) :
        osg::Referenced(),
        N(3),
        triShape(s),
        R0(R0_),
        oneOverR0(1.0/R0) {
        const orsa::TriShape::VertexVector & vv = triShape->getVertexVector();
        const orsa::TriShape::FaceVector   & fv = triShape->getFaceVector();
        aux.resize(fv.size());
        for (size_t fi=0; fi<fv.size(); ++fi) {
            aux[fi].simplexVertexVector.resize(N+1);
            aux[fi].simplexVertexVector[0] = oneOverR0*orsa::Vector(0,0,0);
            aux[fi].simplexVertexVector[1] = oneOverR0*vv[fv[fi].i()];
            aux[fi].simplexVertexVector[2] = oneOverR0*vv[fv[fi].j()];
            aux[fi].simplexVertexVector[3] = oneOverR0*vv[fv[fi].k()];
            //
            // volume of simplex with the face as base and the origin as 4th vertex
            // if moving 4th point away from origin, more terms must be added
            // aux[fi].volume = fabs((vv[fv[fi].i()]*orsa::externalProduct(vv[fv[fi].j()],vv[fv[fi].k()])) / 6);
            //
            // generic, just in case
            aux[fi].volume =
                fabs((aux[fi].simplexVertexVector[1]-aux[fi].simplexVertexVector[0]) *
                     orsa::externalProduct(aux[fi].simplexVertexVector[2]-aux[fi].simplexVertexVector[0],
                                           aux[fi].simplexVertexVector[3]-aux[fi].simplexVertexVector[0]) / 6.0);
            // ORSA_DEBUG("fi: %02i  volume: %g",fi,(*aux[fi].volume));
        }
    }
protected:
    virtual ~SimplexIntegration() { }

public:
    void reserve(const size_t maxDegree) {
        const size_t maxIndex = std::max(getIndex(maxDegree,0,0),
                                         std::max(getIndex(0,maxDegree,0),
                                                  getIndex(0,0,maxDegree)));
        val.reserve(maxIndex+1);
        // val_vol_sum_fun.reserve(maxIndex+1); 
    }
protected:
    class SimplexInternals {
    public:
        std::vector<orsa::Vector> simplexVertexVector;
        orsa::Cache<mpf_class> volume;
    };
protected:
    const size_t N; // N = 3 = dimension of space = number of vertexes in n-dim simplex + 1
    osg::ref_ptr<const orsa::TriShape> triShape;
    const double R0;
    const double oneOverR0;
    mutable std::vector< orsa::Cache<double> > val; // integral value
    mutable std::vector< SimplexInternals > aux;
    // val_vol_sum_fun is the sum over all simplexes of the funciton of given q times the volume of each simplex
    // mutable std::vector< std::vector< orsa::Cache<mpf_class> > > val_vol_sum_fun;
public:
    double getIntegral(const size_t & nx, const size_t & ny, const size_t & nz) const {
        const size_t degree = nx+ny+nz;
        const size_t index = getIndex(nx,ny,nz);
        if (val.size() <= index) {
            val.resize(index+1);
        }
        if (!val[index].isSet()) {
            std::vector<mpf_class> val_vol_sum_fun;
            val_vol_sum_fun.resize(degree+1);
            const size_t q_min = (degree==0) ? 0 : 1;
            const orsa::TriShape::FaceVector & fv = triShape->getFaceVector();
            std::vector<size_t> indexVector;
            for (size_t q=q_min; q<=degree; ++q) {
                mpf_class sum_vol_fun("0.0");
                indexVector.resize(q);
                for (size_t fi=0; fi<fv.size(); ++fi) {
                    for (size_t i=0; i<q; ++i) {
                        indexVector[i] = 0;
                    }
                    mpf_class sum_vol_fi("0.0");
                    while (1) {
                        // orsa::Vector vI(0,0,0);
                        mpf_class vIx = mpf_class("0.0");
                        mpf_class vIy = mpf_class("0.0");
                        mpf_class vIz = mpf_class("0.0");
                        for (size_t i=0; i<q; ++i) {
                            // vI += aux[fi].simplexVertexVector[indexVector[i]];
                            // ORSA_DEBUG("iV[%02i] = %i",i,indexVector[i]);
                            vIx += aux[fi].simplexVertexVector[indexVector[i]].getX();
                            vIy += aux[fi].simplexVertexVector[indexVector[i]].getY();
                            vIz += aux[fi].simplexVertexVector[indexVector[i]].getZ();
                        }
                        /* sum_vol_fi +=
                           orsa::int_pow(vI.getX(),nx)*
                           orsa::int_pow(vI.getY(),ny)*
                           orsa::int_pow(vI.getZ(),nz);
                        */
                        mpf_class vI_pow("1.0");
                        for (size_t p=0; p<nx; ++p) { vI_pow *= vIx; }
                        for (size_t p=0; p<ny; ++p) { vI_pow *= vIy; }
                        for (size_t p=0; p<nz; ++p) { vI_pow *= vIz; }
                        sum_vol_fi += vI_pow;
                        /* ORSA_DEBUG("vI = %+g %+g %+g sum term: %+16.6e",vI.getX(),vI.getY(),vI.getZ(),
                           orsa::int_pow(vI.getX(),nx)*
                           orsa::int_pow(vI.getY(),ny)*
                           orsa::int_pow(vI.getZ(),nz));
                        */
                        bool increased = false;
                        for (size_t i=0; i<q; ++i) {
                            if (indexVector[i]<N) {
                                ++indexVector[i];
                                increased = true;
                                for (size_t s=0; s<i; ++s) {
#warning IMPORTANT: which rule is correct?
                                    indexVector[s] = indexVector[i]; // this one avoids repetitions, i.e. {1,0}, {0,1}
                                    // indexVector[s] = 0; // this one includes repetitions
                                }
                                break;
                            }
                        }
                        // ORSA_DEBUG("increased: %i",increased);
                        if (!increased) {
                            break;
                        }                                
                    }
                    sum_vol_fun += sum_vol_fi*(*aux[fi].volume);
                    /* ORSA_DEBUG("degree: %02i  q: %02i sum_vol_fi: %+16.6e  partial sum_vol_fun: %+16.6e",
                       degree,q,sum_vol_fi,sum_vol_fun);
                    */
                }
                val_vol_sum_fun[q] = sum_vol_fun;
            }
            mpf_class retVal("0.0");
            for (size_t q=q_min; q<=degree; ++q) {
                // const double sum_vol_fun_factor = orsa::binomial(N+degree,N+q).get_d();
                // const int sign = orsa::power_sign(degree-q);
                // retVal += sign*sum_vol_fun_factor*(*val_vol_sum_fun[q]);
                retVal +=
                    orsa::power_sign(degree-q) *
                    orsa::binomial(N+degree,N+q) *
                    val_vol_sum_fun[q];
                /* ORSA_DEBUG("degree: %i  q: %i  factor = binomial(%i,%i): : %g  sign: %i  term: %g",
                   degree,
                   q,
                   N+degree,
                   N+q,
                   sum_vol_fun_factor,sign,
                   val_vol_sum_fun[q]);
                */
            }
            // val[index] = mpf_class(retVal / orsa::binomial(N+degree,degree) / orsa::factorial(degree)).get_d();
            val[index] = mpf_class(retVal / (orsa::binomial(N+degree,degree)*orsa::factorial(degree))).get_d();
        }
        return val[index];
    }
public:
    typedef std::vector< std::vector< std::vector<size_t> > > IndexTableType;
    static IndexTableType indexTable;
public:
    static void updateIndexTable(const size_t & requestedDegree) {
        if (indexTable.size() < (requestedDegree+1)) {
            indexTable.resize(requestedDegree+1);
            for (size_t i=0; i<=requestedDegree; ++i) {
                indexTable[i].resize(requestedDegree+1-i);
                for (size_t j=0; j<=requestedDegree-i; ++j) {
                    indexTable[i][j].resize(requestedDegree+1-i-j);
                }
            }
            size_t idx=0;
            size_t degree=0;
            while (degree <= requestedDegree) {
                for (unsigned int i=0; i<=degree; ++i) {
                    for (unsigned int j=0; j<=degree; ++j) {
                        for (unsigned int k=0; k<=degree; ++k) {
                            if (i+j+k==degree) {
                                // ORSA_DEBUG("inserting %i-%i-%i  index: %i",i,j,k,idx);
                                indexTable[i][j][k] = idx++;
                            }
                        }
                    }
                }
                ++degree;
            }
        }
    }
public:
    static size_t getIndex(const size_t & nx, const size_t & ny, const size_t & nz) {
        const size_t requestedDegree=nx+ny+nz;
        updateIndexTable(requestedDegree);
        return indexTable[nx][ny][nz]; 
    }

    // same as above, but for 4 indices
public:
    typedef std::vector< std::vector< std::vector< std::vector<size_t> > > > Index4TableType;
    static Index4TableType index4Table;
public:
    static void updateIndex4Table(const size_t & requestedDegree) {
        if (index4Table.size() < (requestedDegree+1)) {
            index4Table.resize(requestedDegree+1);
            for (size_t i=0; i<=requestedDegree; ++i) {
                index4Table[i].resize(requestedDegree+1-i);
                for (size_t j=0; j<=requestedDegree-i; ++j) {
                    index4Table[i][j].resize(requestedDegree+1-i-j);
                    for (size_t k=0; k<=requestedDegree-i-j; ++k) {
                        index4Table[i][j][k].resize(requestedDegree+1-i-j-k);
                    }
                }
            }
            size_t idx=0;
            size_t degree=0;
            while (degree <= requestedDegree) {
                for (unsigned int i=0; i<=degree; ++i) {
                    for (unsigned int j=0; j<=degree; ++j) {
                        for (unsigned int k=0; k<=degree; ++k) {
                            for (unsigned int l=0; l<=degree; ++l) {
                                if (i+j+k+l==degree) {
                                    // ORSA_DEBUG("inserting %i-%i-%i-%i  index4: %i",i,j,k,l,idx);
                                    index4Table[i][j][k][l] = idx++;
                                }
                            }
                        }
                    }
                }
                ++degree;
            }
        }
    }
public:
    static size_t getIndex4(const size_t & i, const size_t & j, const size_t & k, const size_t & l) {
        const size_t requestedDegree=i+j+k+l;
        updateIndex4Table(requestedDegree);
        return index4Table[i][j][k][l]; 
    }
};

#endif // _SIMPLEX_H_
