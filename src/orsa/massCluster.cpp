#include <orsa/massCluster.h>

using namespace orsa;

double MassCluster::gravitationalPotential(const orsa::MassCluster * M1,
                                           const orsa::Matrix      & A1_g2l,
                                           const orsa::MassCluster * M2,
                                           const orsa::Matrix      & A2_g2l,
                                           const orsa::Vector      & R) {
    ORSA_DEBUG("code needed here...");
    return 0.0;
}

orsa::Vector MassCluster::gravitationalForce(const orsa::MassCluster * M1,
                                             const orsa::Matrix      & A1_g2l,
                                             const orsa::MassCluster * M2,
                                             const orsa::Matrix      & A2_g2l,
                                             const orsa::Vector      & R) {
    orsa::Matrix tmpMatrix;
    orsa::Matrix::transpose(A1_g2l,tmpMatrix);
    const orsa::Matrix A1_l2g = tmpMatrix;
    orsa::Matrix::transpose(A2_g2l,tmpMatrix);
    const orsa::Matrix A2_l2g = tmpMatrix;
    const orsa::MassCluster::MassClusterVector & mcv1 = M1->massClusterVector;
    const orsa::MassCluster::MassClusterVector & mcv2 = M2->massClusterVector;
    const double eps_sq =
        orsa::square(M1->softeningDistance) +
        orsa::square(M2->softeningDistance);
    orsa::Vector F(0,0,0);
    for (size_t k1=0; k1<mcv1.size(); ++k1) {
        for (size_t k2=0; k2<mcv2.size(); ++k2) {
            const orsa::Vector d12 =
                R + A2_l2g*mcv2[k2].position - A1_l2g*mcv1[k1].position;
            F += mcv2[k2].mass*d12/pow(d12.lengthSquared()+eps_sq,1.5);
        }
    }
    return F;
}

orsa::Vector MassCluster::gravitationalTorque(const orsa::MassCluster * M1,
                                              const orsa::Matrix      & A1_g2l,
                                              const orsa::MassCluster * M2,
                                              const orsa::Matrix      & A2_g2l,
                                              const orsa::Vector      & R) {
    ORSA_DEBUG("code needed here...");
    return orsa::Vector(0,0,0);
}
