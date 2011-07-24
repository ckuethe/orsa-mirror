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
protected:
    mutable std::vector< orsa::Cache<double> > val;
public:
    double getIntegral(const size_t & nx, const size_t & ny, const size_t & nz) const {
        const size_t index = getIndex(nx,ny,nz);
        if (val.size() <= index) {
            val.resize(index+1);
            // compute it
            const orsa::TriShape::VertexVector vv = triShape->getVertexVector();
            const orsa::TriShape::FaceVector   fv = triShape->getFaceVector();
            double sum = 0.0;
            for (size_t fi=0; fi<fv.size(); ++fi) {
                // simplex by simplex...
                // volume of simplex with the face as base and the origin as 4th vertex
                const double volume = (vv[fv[fi].i()]*orsa::externalProduct(vv[fv[fi].j()],vv[fv[fi].k()])) / 6;
                
                // more code here...
                ORSA_DEBUG("more code required here!!");
                
                sum += volume; // / (...) * (...)
                
            }
            val[index] = sum;
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
