#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/shape.h>

class SimplexIntegration : public osg::Referenced {
public:
    SimplexIntegration(const orsa::TriShape * s) :
        osg::Referenced(),
        triShape(s) {
        const orsa::TriShape::VertexVector & vv = triShape->getVertexVector();
        const orsa::TriShape::FaceVector   & fv = triShape->getFaceVector();
        aux.resize(fv.size());
        for (size_t fi=0; fi<fv.size(); ++fi) {
            aux[fi].simplexVertexVector.resize(4);
            aux[fi].simplexVertexVector[0] = orsa::Vector(0,0,0);
            aux[fi].simplexVertexVector[1] = vv[fv[fi].i()];
            aux[fi].simplexVertexVector[2] = vv[fv[fi].j()];
            aux[fi].simplexVertexVector[3] = vv[fv[fi].k()];
            //
            // volume of simplex with the face as base and the origin as 4th vertex
            // if moving 4th point away from origin, more terms must be added
            // aux[fi].volume = (vv[fv[fi].i()]*orsa::externalProduct(vv[fv[fi].j()],vv[fv[fi].k()])) / 6;
            // generic, just in case
            aux[fi].volume =
                (aux[fi].simplexVertexVector[1]-aux[fi].simplexVertexVector[0]) *
                orsa::externalProduct(aux[fi].simplexVertexVector[2]-aux[fi].simplexVertexVector[0],
                                      aux[fi].simplexVertexVector[3]-aux[fi].simplexVertexVector[0]) / 6.0;
        }
    }
protected:
    virtual ~SimplexIntegration() { }
protected:
    osg::ref_ptr<const orsa::TriShape> triShape;
protected:
    class SimplexInternals {
    public:
        // val_p[index][index4]...
        mutable std::vector< std::vector< orsa::Cache<double> > > val_p;
        std::vector<orsa::Vector> simplexVertexVector;
        orsa::Cache<double> volume;
    };
public:
    static double H(const size_t & nx, const size_t & ny, const size_t & nz,
                    const std::vector<size_t> indexVector,
                    SimplexInternals & aux) {
        const size_t degree = nx+ny+nz;
        // ORSA_DEBUG("degree: %d",degree);
        std::vector<size_t> i4;
        i4.resize(4);
        for (size_t s=0; s<4; ++s) {
            i4[s] = 0;
        }
        for (size_t q=0; q<degree; ++q) {
            ++i4[indexVector[q]];
        }
        /* for (size_t s=0; s<4; ++s) {
           ORSA_DEBUG("i4[%d] = %d",s,i4[s]);
           }
        */
        const size_t index  = getIndex(nx,ny,nz);
        const size_t index4 = getIndex4(i4[0],i4[1],i4[2],i4[3]);
        if (aux.val_p.size() <= index) {
            aux.val_p.resize(index+1);
            if (aux.val_p[index].size() <= index4) {
                aux.val_p[index].resize(index4+1);
            }
        }
        // ORSA_DEBUG("val_p[%i][%i] set: %i",index,index4,aux.val_p[index][index4].isSet());
        const size_t maxCount = 1<<degree; // 2^degree
        std::vector<bool> I;
        I.resize(degree);
        std::vector<orsa::Vector> localSimplexVertexVector;
        localSimplexVertexVector.resize(degree);
        for (size_t q=0; q<degree; ++q) {
            localSimplexVertexVector[q] = aux.simplexVertexVector[indexVector[q]];
        }
        std::vector<size_t> t4;
        t4.resize(4);
        double retVal = 0.0;
        for (size_t count=0; count<maxCount; ++count) {
            for (size_t s=0; s<4; ++s) {
                t4[s] = 0;
            }
            for (size_t i=0; i<degree; ++i) {
                I[i] = (count & (1<<i));
                if (I[i]) {
                    ++t4[indexVector[i]];
                }
                // ORSA_DEBUG("count: %d   I[%d] = %i",count,i,I[i]==true);
            }
            /* for (size_t s=0; s<4; ++s) {
               ORSA_DEBUG("t4[%d] = %d",s,t4[s]);
               }
            */
            const size_t it4 = getIndex4(t4[0],t4[1],t4[2],t4[3]);
            if (aux.val_p[index].size() <= it4) {
                // is this needed?
                aux.val_p[index].resize(it4+1);
            }
            const int sign = orsa::power_sign(degree-(t4[0]+t4[1]+t4[2]+t4[3]));
            if (aux.val_p[index][it4].isSet()) {
                retVal += aux.val_p[index][it4];
                // ORSA_DEBUG("--cached--");
            } else {
                // ORSA_DEBUG("--NOT-cached--");
                orsa::Vector vertex_I(0,0,0);
                size_t num_I = 0;
                for (size_t i=0; i<degree; ++i) {
                    if (I[i]) {
                        vertex_I += localSimplexVertexVector[i];  
                        ++num_I;
                    }
                }
                // const int sign = orsa::power_sign(degree-num_I);
                const double p =
                    orsa::int_pow(vertex_I.getX(),nx)*
                    orsa::int_pow(vertex_I.getY(),ny)*
                    orsa::int_pow(vertex_I.getZ(),nz);
                aux.val_p[index][it4] = p;
                retVal += sign*p;
                // ORSA_DEBUG("sign: %+i   p: %16.6e",sign,p);
            }
        }
        retVal /= orsa::factorial(nx+ny+nz).get_d();
        // ORSA_DEBUG("retVal: %16.6e",retVal);
        // aux.val_p[index][index4] = retVal;
        // ORSA_DEBUG("setting val_p[%d][%d]",index,index4);
        return retVal;                
        // return aux.val_p[index][index4];
        // return retVal;
    }
public:
    static double sum_H(const size_t & nx, const size_t & ny, const size_t & nz,
                        SimplexInternals & aux) {
        const size_t N = aux.simplexVertexVector.size(); // N here is 4
        const size_t degree = nx+ny+nz;
        
        std::vector<size_t> indexVector;
        indexVector.resize(degree);
        
        for (size_t q=0; q<degree; ++q) {
            indexVector[q] = 0;
        }
        
        // std::vector<orsa::Vector> simplexVertexVector_Degree;
        // simplexVertexVector_Degree.resize(degree);
        
#warning maybe use mpf_class for retVal?
        double retVal = 0.0;
        // size_t calls = 0;
        bool done = false;
        do {
            // first iter with indexVector={0,0...,0}
            
            /* for (size_t q=0; q<degree; ++q) {
               ORSA_DEBUG("indexVector[%d] = %i",q,indexVector[q]);
               }
            */
            
            /* for (size_t q=0; q<degree; ++q) {
               simplexVertexVector_Degree[q] = simplexVertexVector[indexVector[q]];
               }
            */
            
            // retVal += H(nx,ny,nz,simplexVertexVector_Degree);    
            retVal += H(nx,ny,nz,indexVector,aux);
            
            // ++calls;
            
            bool increased = false;
            for (size_t q=0; q<degree; ++q) {
                if (indexVector[q]<(N-1)) {
                    ++indexVector[q];
                    increased = true;
                    for (size_t s=0; s<q; ++s) {
                        indexVector[s] = indexVector[q];
                    }
                    break;
                }
            }
            // ORSA_DEBUG("increased: %i",increased);
            
            if (!increased) {
                done = true;
            }
            
        } while (!done);
        
        return retVal;
    }
protected:
    mutable std::vector< orsa::Cache<double> > val;
    mutable std::vector< SimplexInternals > aux;
public:
    double getIntegral(const size_t & nx, const size_t & ny, const size_t & nz) const {
        const size_t degree = nx+ny+nz;
        const size_t index = getIndex(nx,ny,nz);
        if (val.size() <= index) {
            val.resize(index+1);
        }
        if (!val[index].isSet()) {
            // compute it
            // const orsa::TriShape::VertexVector & vv = triShape->getVertexVector();
            const orsa::TriShape::FaceVector   & fv = triShape->getFaceVector();
            // std::vector<orsa::Vector> simplexVertexVector;
            // simplexVertexVector.resize(4);
            double sum = 0.0;
            for (size_t fi=0; fi<fv.size(); ++fi) {
                // simplex by simplex...
                
#warning if shape is strongly concave and a simplex covers volume outside the body shape, then the results are incorrect (including volume computations...)
                
#warning default origin for 4th simplex vertex, should be a parameter of the class??                
                
                sum += aux[fi].volume*sum_H(nx,ny,nz,aux[fi]);
            }
            val[index] = sum / orsa::binomial(3+degree,degree).get_d();
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
