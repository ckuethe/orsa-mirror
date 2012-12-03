#ifndef _ORSA_MASS_CLUSTER_H_
#define _ORSA_MASS_CLUSTER_H_

#include <vector>

#include <orsa/massDistribution.h>
#include <orsa/matrix.h>
#include <orsa/shape.h>
#include <orsa/vector.h>

#include <osg/Referenced>
#include <osg/ref_ptr>

namespace orsa {
    
    class MassCluster : public osg::Referenced {
        class MassClusterElement {
        public:
            orsa::Vector position;
            double mass;
        };
        typedef std::vector<MassClusterElement> MassClusterVector;
    public:
        MassCluster(const orsa::Shape * shape,
                    const orsa::MassDistribution * massDistribution,
                    const size_t & samplePoints,
                    const double & softeningDistance_,
                    const orsa::Cache<orsa::Vector> nominalCenterOfMass);
    protected:
        virtual ~MassCluster() { }
    public:
        const double softeningDistance;
        MassClusterVector massClusterVector;
    public:
        static double gravitationalPotential(const orsa::MassCluster * M1,
                                             const orsa::Matrix      & A1_g2l,
                                             const orsa::MassCluster * M2,
                                             const orsa::Matrix      & A2_g2l,
                                             const orsa::Vector      & R);
    public:
        static orsa::Vector gravitationalForce(const orsa::MassCluster * M1,
                                               const orsa::Matrix      & A1_g2l,
                                               const orsa::MassCluster * M2,
                                               const orsa::Matrix      & A2_g2l,
                                               const orsa::Vector      & R);
    public:
        static orsa::Vector gravitationalTorque(const orsa::MassCluster * M1,
                                                const orsa::Matrix      & A1_g2l,
                                                const orsa::MassCluster * M2,
                                                const orsa::Matrix      & A2_g2l,
                                                const orsa::Vector      & R);
    };
    
}; // namespace orsa

#endif // _ORSA_MASS_CLUSTER_H_
