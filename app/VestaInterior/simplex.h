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
        
    }
protected:
    virtual ~SimplexIntegration() { }
protected:
    osg::ref_ptr<const orsa::TriShape> triShape;
public:
    static double H(const size_t & nx, const size_t & ny, const size_t & nz,
                    const std::vector<orsa::Vector> & simplexVertexVector) {
        const size_t N = simplexVertexVector.size();
        const size_t maxCount = 1<<N; // 2^N
        std::vector<bool> I;
        I.resize(N);
        double retVal = 0.0;
        for (size_t count=0; count<maxCount; ++count) {
            for (size_t i=0; i<N; ++i) {
                I[i] = (count & (1<<i));
                // ORSA_DEBUG("count: %d   I[%d] = %i",count,i,I[i]==true);
            }
            orsa::Vector vertex_I(0,0,0);
            size_t num_I = 0;
            for (size_t i=0; i<N; ++i) {
                if (I[i]) {
                    vertex_I += simplexVertexVector[i];  
                    ++num_I;
                }
            }
            const int sign = orsa::power_sign(N-num_I);
            const double p =
                orsa::int_pow(vertex_I.getX(),nx)*
                orsa::int_pow(vertex_I.getY(),ny)*
                orsa::int_pow(vertex_I.getZ(),nz);
            retVal += sign*p;
            // ORSA_DEBUG("sign: %+i   p: %16.6e",sign,p);
        }
        retVal /= orsa::factorial(nx+ny+nz).get_d();
        // ORSA_DEBUG("retVal: %16.6e",retVal);
        return retVal;
    }
public:
    static double sum_H(const size_t & nx, const size_t & ny, const size_t & nz,
                        const std::vector<orsa::Vector> & simplexVertexVector) {
        const size_t N = simplexVertexVector.size();
        const size_t degree = nx+ny+nz;

        mpz_class z = 1;
        for (size_t q=0; q<degree; ++q) {
            z *= N;
        }
        const mpz_class maxCount = z; // N^q

        std::vector<size_t> indexVector;
        indexVector.resize(degree);
        
        std::vector<orsa::Vector> simplexVertexVector_Degree;
        simplexVertexVector_Degree.resize(degree);
        
#warning maybe use mpf_class? (at least internally, before returning value?)
        double retVal = 0.0;
        for (mpz_class z=0; z<maxCount; ++z) {
            for (size_t q=0; q<degree; ++q) {
                mpz_class local_z = z;
                for (size_t s=0; s<q; ++s) {
                    local_z /= N;
                }
                indexVector[q] = mpz_class(local_z % N).get_ui();
                // ORSA_DEBUG("z: %Zd   indexVector[%d] = %i",z.get_mpz_t(),q,indexVector[q]);
            }
            
            for (size_t q=0; q<degree; ++q) {
                simplexVertexVector_Degree[q] = simplexVertexVector[indexVector[q]];
            }
            
            retVal += H(nx,ny,nz,simplexVertexVector_Degree);           
        }
        return retVal;
    }
protected:
    mutable std::vector< orsa::Cache<double> > val;
public:
    double getIntegral(const size_t & nx, const size_t & ny, const size_t & nz) const {
        const size_t degree = nx+ny+nz;
        const size_t index = getIndex(nx,ny,nz);
        if (val.size() <= index) {
            val.resize(index+1);
        }
        if (!val[index].isSet()) {
            // compute it
            const orsa::TriShape::VertexVector & vv = triShape->getVertexVector();
            const orsa::TriShape::FaceVector   & fv = triShape->getFaceVector();
            std::vector<orsa::Vector> simplexVertexVector;
            simplexVertexVector.resize(4);
            double sum = 0.0;
            for (size_t fi=0; fi<fv.size(); ++fi) {
                // simplex by simplex...
                // volume of simplex with the face as base and the origin as 4th vertex
                // if moving 4th point away from origin, more terms must be added
                // also if shape is strongly concave and a simplex covers volume outside the body shape, then the results are incorrect
                const double volume = (vv[fv[fi].i()]*orsa::externalProduct(vv[fv[fi].j()],vv[fv[fi].k()])) / 6;
                
#warning ONE SUMMATION IS MISSING HERE, either on sum alpha_i = q OR on 0 <= i1 <= i2 ... <= iq <= n
#warning ALSO, MUST RESIZE simplexVertexVector, from (dim+1)=4 to the degree!! (q-homogeneous...)
                
#warning default origin for 4th simplex vertex, should be a parameter of the class??
                
                
#warning RESIZE to DEGREE, and do all the combinations!
                simplexVertexVector[0] = orsa::Vector(0,0,0);
                simplexVertexVector[1] = vv[fv[fi].i()];
                simplexVertexVector[2] = vv[fv[fi].j()];
                simplexVertexVector[3] = vv[fv[fi].k()];
                
                
                // sum += volume*H(nx,ny,nz,simplexVertexVector);
                sum += volume*sum_H(nx,ny,nz,simplexVertexVector);
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
};

#endif // _SIMPLEX_H_
