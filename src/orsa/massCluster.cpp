#include <orsa/massCluster.h>

#include <orsa/util.h>

using namespace orsa;

MassCluster::MassCluster(const orsa::Shape * shape,
                         const orsa::MassDistribution * massDistribution,
                         const size_t & samplePoints,
                         const orsa::Cache<orsa::Vector> nominalCenterOfMass) :
    osg::Referenced(true) {
    
    if (!shape) return;
    if (!massDistribution) return;
    
    osg::ref_ptr<RandomPointsInShape> randomPointsInShape =
        new RandomPointsInShape(shape,
                                massDistribution,
                                samplePoints,
                                false);
    const orsa::Vector centerOfMass = orsa::centerOfMass(randomPointsInShape);
    orsa::Vector global_offset(0,0,0);
    if (nominalCenterOfMass.isSet()) {
        global_offset = nominalCenterOfMass - centerOfMass;
        // ORSA_DEBUG("global offset: %g [km]",orsa::FromUnits(global_offset.length(),orsa::Unit::KM,-1));
    }
    orsa::Vector v;
    double density;
    massClusterVector.reserve(samplePoints);
    randomPointsInShape->reset();
    while (randomPointsInShape->get(v,density)) {
        if (density > 0) {
            // enforce center of mass position?
            if (nominalCenterOfMass.isSet()) {
                // v += nominalCenterOfMass - centerOfMass;
                v += global_offset;
            }
            MassCluster::MassClusterElement el;
            el.position = v;
            el.density  = density;
            massClusterVector.push_back(el);
        }
    }
    if (massClusterVector.size() == 0) {
        ORSA_DEBUG("problems...");
        exit(0);
    }
    // softeningDistance = cbrt(shape->volume()/massClusterVector.size());
    // ORSA_DEBUG("softeningDistance: %g [km]",orsa::FromUnits(softeningDistance,orsa::Unit::KM,-1));
    //
    vR = cbrt(shape->volume()/(massClusterVector.size()*(4.0*orsa::pi()/3.0)));
    vR2 = vR*vR;
    vR3 = vR2*vR;
    // ORSA_DEBUG("vR: %g [km]",orsa::FromUnits(vR,orsa::Unit::KM,-1));
}

void MassCluster::print() const {
    ORSA_DEBUG("vR: %g [km]",orsa::FromUnits(vR,orsa::Unit::KM,-1));
    ORSA_DEBUG("massClusterVector.size(): %i",massClusterVector.size());
    for (size_t k=0; k<massClusterVector.size(); ++k) {
        const MassCluster::MassClusterElement el = massClusterVector[k];
        ORSA_DEBUG("element[%06i] density: %g [g/cm^3]   position: %g %g %g [km]",
                   k,
                   orsa::FromUnits(orsa::FromUnits(el.density,orsa::Unit::GRAM,-1),orsa::Unit::CM,3),
                   orsa::FromUnits(el.position.getX(),orsa::Unit::KM,-1),
                   orsa::FromUnits(el.position.getY(),orsa::Unit::KM,-1),
                   orsa::FromUnits(el.position.getZ(),orsa::Unit::KM,-1));
    }
}

double MassCluster::gravitationalPotential(const orsa::MassCluster * M1,
                                           const orsa::Matrix      & A1_g2l,
                                           const orsa::MassCluster * M2,
                                           const orsa::Matrix      & A2_g2l,
                                           const orsa::Vector      & R) {
                                               ORSA_DEBUG("code needed...");
                                               return 0.0;
                                               /*
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
    double V=0.0;
    double sum_weight=0.0;
    for (size_t k1=0; k1<mcv1.size(); ++k1) {
        for (size_t k2=0; k2<mcv2.size(); ++k2) {
            const orsa::Vector d12 =
                R + A2_l2g*mcv2[k2].position - A1_l2g*mcv1[k1].position;
            const double weight = mcv1[k1].density * mcv2[k2].density;
            V -= weight/pow(d12.lengthSquared()+eps_sq,0.5);
            sum_weight += weight;
        }
    }
    if (sum_weight != 0.0) V /= sum_weight;
    // ORSA_DEBUG("R: %g   1/R^2 = %g   F = %g",R.length(),1.0/R.lengthSquared(),F.length());
    // ORSA_DEBUG("QQQ R: %g acos(R*F): %g",R.length(),orsa::radToDeg()*acos(R.normalized()*F.normalized()));
    return V;    
                                               */                                       
}

double MassCluster::gravitationalPotential(const orsa::MassCluster * M1,
                                           const orsa::Matrix      & A1_g2l,
                                           // const orsa::MassCluster * M2,
                                           // const orsa::Matrix      & A2_g2l,
                                           const orsa::Vector      & R) {
    orsa::Matrix tmpMatrix;
    orsa::Matrix::transpose(A1_g2l,tmpMatrix);
    const orsa::Matrix A1_l2g = tmpMatrix;
    // orsa::Matrix::transpose(A2_g2l,tmpMatrix);
    // const orsa::Matrix A2_l2g = tmpMatrix;
    const orsa::MassCluster::MassClusterVector & mcv1 = M1->massClusterVector;
    // const orsa::MassCluster::MassClusterVector & mcv2 = M2->massClusterVector;
    // const double eps_sq = orsa::square(M1->softeningDistance); // +
        // orsa::square(M2->softeningDistance);
    const double & vR  = M1->vR;
    const double & vR2 = M1->vR2;
    const double & vR3 = M1->vR3;
    double V=0.0;
    double sum_weight=0.0;
    for (size_t k1=0; k1<mcv1.size(); ++k1) {
        // for (size_t k2=0; k2<mcv2.size(); ++k2) {
            const orsa::Vector d12 = R - A1_l2g*mcv1[k1].position;
            // R + A2_l2g*mcv2[k2].position - A1_l2g*mcv1[k1].position;
            const double weight = mcv1[k1].density;
            // const double weight = mcv1[k1].density * mcv2[k2].density;
            // V -= weight/pow(d12.lengthSquared()+eps_sq,0.5);
            // V -= weight/sqrt(d12.lengthSquared()+eps_sq);
            if (d12.lengthSquared()>vR2) {
                V -= weight/d12.length();
            } else {
                V -= weight*(3*vR2-d12.lengthSquared())/(2*vR3);
            }
            sum_weight += weight;
        // }
    }
    if (sum_weight != 0.0) V /= sum_weight;
    // ORSA_DEBUG("R: %g   1/R^2 = %g   F = %g",R.length(),1.0/R.lengthSquared(),F.length());
    // ORSA_DEBUG("QQQ R: %g acos(R*F): %g",R.length(),orsa::radToDeg()*acos(R.normalized()*F.normalized()));
    return V;                                           
}

orsa::Vector MassCluster::gravitationalForce(const orsa::MassCluster * M1,
                                             const orsa::Matrix      & A1_g2l,
                                             const orsa::MassCluster * M2,
                                             const orsa::Matrix      & A2_g2l,
                                             const orsa::Vector      & R) {
                                                 ORSA_DEBUG("code needed...");
                                                 return orsa::Vector(0,0,0);
                                                 /*
                                                 
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
    double sum_weight=0.0;
    for (size_t k1=0; k1<mcv1.size(); ++k1) {
        for (size_t k2=0; k2<mcv2.size(); ++k2) {
            const orsa::Vector d12 =
                R + A2_l2g*mcv2[k2].position - A1_l2g*mcv1[k1].position;
            const double weight = mcv1[k1].density * mcv2[k2].density;
            F += weight*d12/pow(d12.lengthSquared()+eps_sq,1.5);
            sum_weight += weight;
        }
    }
    if (sum_weight != 0.0) F /= sum_weight;
    // ORSA_DEBUG("R: %g   1/R^2 = %g   F = %g",R.length(),1.0/R.lengthSquared(),F.length());
    // ORSA_DEBUG("QQQ R: %g acos(R*F): %g",R.length(),orsa::radToDeg()*acos(R.normalized()*F.normalized()));
    return F;
                                                 */
}

orsa::Vector MassCluster::gravitationalAcceleration(const orsa::MassCluster * M1,
                                             		const orsa::Matrix      & A1_g2l,
                                             	   	// const orsa::MassCluster * M2,
                                             	   	// const orsa::Matrix      & A2_g2l,
                                             	   	const orsa::Vector      & R) {
    orsa::Matrix tmpMatrix;
    orsa::Matrix::transpose(A1_g2l,tmpMatrix);
    const orsa::Matrix A1_l2g = tmpMatrix;
    // orsa::Matrix::transpose(A2_g2l,tmpMatrix);
    // const orsa::Matrix A2_l2g = tmpMatrix;
    const orsa::MassCluster::MassClusterVector & mcv1 = M1->massClusterVector;
    // const orsa::MassCluster::MassClusterVector & mcv2 = M2->massClusterVector;
    // const double eps_sq = orsa::square(M1->softeningDistance); // + // orsa::square(M2->softeningDistance);
    const double & vR  = M1->vR;
    const double & vR2 = M1->vR2;
    const double & vR3 = M1->vR3;
    orsa::Vector F(0,0,0);
    double sum_weight=0.0;
    for (size_t k1=0; k1<mcv1.size(); ++k1) {
        // for (size_t k2=0; k2<mcv2.size(); ++k2) {
            const orsa::Vector d12 = R - A1_l2g*mcv1[k1].position;
            // R + A2_l2g*mcv2[k2].position - A1_l2g*mcv1[k1].position;
        	const double weight = mcv1[k1].density;
        	// const double weight = mcv1[k1].density * mcv2[k2].density;
            // F += weight*d12/pow(d12.lengthSquared()+eps_sq,1.5);
            // F += weight*d12/pow(d12.lengthSquared()+eps_sq,1.5);
            if (d12.lengthSquared()>vR2) {
                F += weight*d12/orsa::cube(d12.length());
            } else {
                F += weight*d12/vR3;
            }
            sum_weight += weight;
		// }
    }
    if (sum_weight != 0.0) F /= sum_weight;
    // ORSA_DEBUG("R: %g   1/R^2 = %g   F = %g",R.length(),1.0/R.lengthSquared(),F.length());
    // ORSA_DEBUG("QQQ R: %g acos(R*F): %g",R.length(),orsa::radToDeg()*acos(R.normalized()*F.normalized()));
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
