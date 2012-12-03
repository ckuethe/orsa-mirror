#ifndef _ORSA_MASS_CLUSTER_H_
#define _ORSA_MASS_CLUSTER_H_

#include <vector>

#include <orsa/matrix.h>
#include <orsa/vector.h>

#include <osg/Referenced>
#include <osg/ref_ptr>

namespace orsa {
    
    class MassCluster : public osg::Referenced {
        class MassClusterElement {
            orsa::Vector position;
            double mass;
        };
        typedef std::vector<MassClusterElement> MassClusterVector;
    public:
        MassCluster() : osg::Referenced(true) {
            softeningDistance=0.0;
        }
    protected:
        virtual ~MassCluster() { }
    public:
        double softeningDistance;
        MassClusterVector massClusterVector;
    };

    static double gravitationalPotential(const orsa::MassCluster * M1,
                                         const orsa::Matrix      & A1_g2l,
                                         const orsa::MassCluster * M2,
                                         const orsa::Matrix      & A2_g2l,
                                         const orsa::Vector      & R);
    
    static orsa::Vector gravitationalForce(const orsa::MassCluster * M1,
                                           const orsa::Matrix      & A1_g2l,
                                           const orsa::MassCluster * M2,
                                           const orsa::Matrix      & A2_g2l,
                                           const orsa::Vector      & R);
    
    static orsa::Vector gravitationalTorque(const orsa::MassCluster * M1,
                                            const orsa::Matrix      & A1_g2l,
                                            const orsa::MassCluster * M2,
                                            const orsa::Matrix      & A2_g2l,
                                            const orsa::Vector      & R);
    
}; // namespace orsa

#endif // _ORSA_MASS_CLUSTER_H_
