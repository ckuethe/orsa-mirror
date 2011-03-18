#include "InsideVesta.h"

#include "vesta.h"

using namespace orsa;
using namespace orsaSolarSystem;

GlobalRNG * GlobalRNG::_instance = 0;

// int main(int argc, char **argv) {
int main() {
    
    orsa::Debug::instance()->initTimer();
    
    // units
    const double    km = orsa::FromUnits(1,orsa::Unit::KM);
    const double g_cm3 = orsa::FromUnits(orsa::FromUnits(1,orsa::Unit::GRAM),orsa::Unit::CM,-3);
    
    const unsigned int degree = 4;
    
    const double vestaMass = 1.35e-10*orsaSolarSystem::Data::MSun();
    /* ORSA_DEBUG("vestaMass = %g x MSun = %.12e kg",
       vestaMass / orsaSolarSystem::Data::MSun(),
       orsa::FromUnits(vestaMass,orsa::Unit::KG,-1));
    */
    
    // change if changing shape
    double volume = FromUnits(7.875e7,Unit::KM,3);
    
    const double bulkDensity = vestaMass / volume;
    ORSA_DEBUG("bulk density: %g [g/cm^3]",
               orsa::FromUnits(orsa::FromUnits(bulkDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
    
    Model model;
    model.totalMass   = new Par(vestaMass,vestaMass);
    model.totalVolume = new Par(volume,volume);
    model.coreDensity = new Par(bulkDensity, 10.0*g_cm3);
    model.coreCenterX = new Par(-50.0*km,    50.0*km);
    model.coreCenterY = new Par(-50.0*km,    50.0*km);
    model.coreCenterZ = new Par(-50.0*km,    50.0*km);
    model.coreRadiusX = new Par(  0.0*km,   200.0*km);
    model.coreRadiusY = new Par(  0.0*km,   200.0*km);
    model.coreRadiusZ = new Par(  0.0*km,   200.0*km);
    
    osg::ref_ptr<orsa::Shape> shape; 
    {
        osg::ref_ptr<VestaShape> vestaShapeThomas = new VestaShape;
        if (!vestaShapeThomas->read("vesta_thomas.dat")) {
            ORSA_ERROR("problems encountered while reading shape file...");
        }
        shape = vestaShapeThomas.get();
    }
    
    // choose here if sample the reference solution...
    // const Model::Values refVal = model.sample();
    // ... or fix it manually (make sure it's within the model limits above!)
    Model::Values tmpVal;
    tmpVal.totalMass   = vestaMass;
    tmpVal.totalVolume = volume;
    tmpVal.coreDensity =   6.5*g_cm3;
    tmpVal.coreCenterX =   2.0*km;
    tmpVal.coreCenterY =   1.0*km;
    tmpVal.coreCenterZ =   0.8*km;
    tmpVal.coreRadiusX = 120.0*km;
    tmpVal.coreRadiusY = 110.0*km;
    tmpVal.coreRadiusZ =  90.0*km;
    const Model::Values refVal = tmpVal;
    
    osg::ref_ptr<ModelMassDistribution> refMD = new ModelMassDistribution(refVal);
    
    if (0) {
        // write profile density, for y=0 plane
        ORSA_DEBUG("writing profile.dat file...");
        const Box boundingBox = shape->boundingBox();
        FILE * fp = fopen("profile.dat","w");
        for (unsigned int k=0; k<1000000; ++k) {
            const orsa::Vector v =  Vector(boundingBox.getXMin()+(boundingBox.getXMax()-boundingBox.getXMin())*GlobalRNG::instance()->gsl_rng_uniform(),
                                           0.0,
                                           boundingBox.getZMin()+(boundingBox.getZMax()-boundingBox.getZMin())*GlobalRNG::instance()->gsl_rng_uniform());
            if (shape->isInside(v)) {
                fprintf(fp,"%g %g %g\n",v.getX(),v.getZ(),refMD->density(v));
            } else {
                fprintf(fp,"%g %g %g\n",v.getX(),v.getZ(),0.0);
            }
        }
        fclose(fp);
        ORSA_DEBUG("done writing profile.dat");
    }
    
    // double dummy_volume;
    orsa::Vector centerOfMass;
    orsa::Matrix shapeToLocal;
    orsa::Matrix localToShape;
    orsa::Matrix inertiaMatrix;
    osg::ref_ptr<orsa::PaulMoment> paulMoment;
    
    const unsigned int N = 1000000;
    // remember, there is another RNG, see GlobalRNG in InsideVesta.h
    // int randomSeed = getpid();
    int randomSeed = 85719;
    
    // save this for later!
    osg::ref_ptr<RandomPointsInShape> randomPointsInShape = new RandomPointsInShape(shape.get(),N,randomSeed);
    
    ORSA_DEBUG("volume from shape: %g   input value: %g   (ratio: %g)",
               orsa::volume(randomPointsInShape),
               volume,
               orsa::volume(randomPointsInShape)/volume);
    
    centerOfMass = orsa::centerOfMass(randomPointsInShape,
                                      refMD.get());
    
    orsa::diagonalizedInertiaMatrix(shapeToLocal,
                                    localToShape,
                                    inertiaMatrix,
                                    centerOfMass,
                                    randomPointsInShape,
                                    refMD.get());
    
    paulMoment = orsa::computePaulMoment(degree,
                                         shapeToLocal,
                                         localToShape,
                                         centerOfMass,
                                         randomPointsInShape,
                                         refMD.get());
    
    const orsa::Vector ref_centerOfMass  = centerOfMass;
    const orsa::Matrix ref_shapeToLocal  = shapeToLocal;
    const orsa::Matrix ref_localToShape  = localToShape;
    const orsa::Matrix ref_inertiaMatrix = inertiaMatrix;
    const osg::ref_ptr<orsa::PaulMoment> ref_paulMoment = paulMoment;
    
    // additional, to track rotational axis offset
    const orsa::Vector ref_rotAxis = ref_localToShape*orsa::Vector(0,0,1);
    
    std::vector< std::vector<double> > ref_C, ref_S, ref_norm_C, ref_norm_S;
    std::vector<double> ref_J;
    orsa::convert(ref_C, ref_S, ref_norm_C, ref_norm_S, ref_J,
                  paulMoment.get(), 
                  FromUnits(300,orsa::Unit::KM));
    
    if (1) {
        // print out...
        std::vector< std::vector<double> > C, S, norm_C, norm_S;
        std::vector<double> J;
        orsa::convert(C, S, norm_C, norm_S, J,
                      paulMoment.get(), 
                      FromUnits(300,orsa::Unit::KM));
        
        ORSA_DEBUG("$\\rho_{m}$ & $%9.3f$ \\\\",orsa::FromUnits(orsa::FromUnits(refMD->_mantleDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
        ORSA_DEBUG("$\\rho_{c}$ & $%9.3f$ \\\\",orsa::FromUnits(orsa::FromUnits(refMD->_val.coreDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
        ORSA_DEBUG("$\\V_{c}$   & $%9.3e$ \\\\",orsa::FromUnits(refMD->_coreVolume,orsa::Unit::KM,-3));
        ORSA_DEBUG("$\\M_{c}$   & $%9.3e$ \\\\",orsa::FromUnits(refMD->_coreVolume*refMD->_val.coreDensity,orsa::Unit::KG,-1));
        // ORSA_DEBUG("$R_{c}$    & $%9.3f$ \\\\", orsa::FromUnits(coreRadius,orsa::Unit::KM,-1));
        ORSA_DEBUG("\%\\hline");
        ORSA_DEBUG("$x_{c}$    & $%+9.3f$ \\\\",orsa::FromUnits(refMD->_val.coreCenterX,orsa::Unit::KM,-1));
        ORSA_DEBUG("$y_{c}$    & $%+9.3f$ \\\\",orsa::FromUnits(refMD->_val.coreCenterY,orsa::Unit::KM,-1));
        ORSA_DEBUG("$z_{c}$    & $%+9.3f$ \\\\",orsa::FromUnits(refMD->_val.coreCenterZ,orsa::Unit::KM,-1));
        ORSA_DEBUG("\%\\hline");
        ORSA_DEBUG("$R_{x}$    & $%+9.3f$ \\\\",orsa::FromUnits(refMD->_val.coreRadiusX,orsa::Unit::KM,-1));
        ORSA_DEBUG("$R_{y}$    & $%+9.3f$ \\\\",orsa::FromUnits(refMD->_val.coreRadiusY,orsa::Unit::KM,-1));
        ORSA_DEBUG("$R_{z}$    & $%+9.3f$ \\\\",orsa::FromUnits(refMD->_val.coreRadiusZ,orsa::Unit::KM,-1));
        ORSA_DEBUG("\\hline");
        ORSA_DEBUG("$x_{0}$    & $%+9.3f$ \\\\",orsa::FromUnits(centerOfMass.getX(),orsa::Unit::KM,-1));
        ORSA_DEBUG("$y_{0}$    & $%+9.3f$ \\\\",orsa::FromUnits(centerOfMass.getY(),orsa::Unit::KM,-1));
        ORSA_DEBUG("$z_{0}$    & $%+9.3f$ \\\\",orsa::FromUnits(centerOfMass.getZ(),orsa::Unit::KM,-1));
        ORSA_DEBUG("\%\\hline");
        for (unsigned int l=2; l<=degree; ++l) {
            // J_l is minus C_l0, where C_l0 is not normalized
            ORSA_DEBUG("$J_{%i}$    & $%+9.6f$ \\\\",l,-C[l][0]);
        }
        ORSA_DEBUG("\%\\hline");
        for (unsigned int l=2; l<=degree; ++l) {
            for (unsigned int m=0; m<=l; ++m) {
                // LaTeX Tabular style
                ORSA_DEBUG("$C_{%i%i}$   & $%+9.6f$ \\\\",l,m,norm_C[l][m]);
                if (m!=0) {
                    ORSA_DEBUG("$S_{%i%i}$   & $%+9.6f$ \\\\",l,m,norm_S[l][m]);
                }
            }
        }
        for (unsigned int l=2; l<=degree; ++l) {
            double var_l=0; // variance_l, or sigma_l squared
            for (unsigned int m=0; m<=l; ++m) {
                var_l += orsa::square(norm_C[l][m])+orsa::square(norm_S[l][m]);
            }
            var_l /= (2*l+1); // definitions may vary...
            ORSA_DEBUG("power: %2i %e",l,var_l);
        }
    }
    
    while (1) {
        const Model::Values val = model.sample();
        osg::ref_ptr<ModelMassDistribution> md = new ModelMassDistribution(val);
        {
            // some consistency checks 
            if (md->_mantleDensity <= 0.0) {
                continue;
            }
        }
        centerOfMass = orsa::centerOfMass(randomPointsInShape,
                                          md.get());
        orsa::diagonalizedInertiaMatrix(shapeToLocal,
                                        localToShape,
                                        inertiaMatrix,
                                        centerOfMass,
                                        randomPointsInShape,
                                        md.get());
        paulMoment = orsa::computePaulMoment(degree,
                                             shapeToLocal,
                                             localToShape,
                                             centerOfMass,
                                             randomPointsInShape,
                                             md.get());
        // derived
        const double centerOfMassOffset = (ref_centerOfMass-centerOfMass).length();
        const orsa::Vector rotAxis = localToShape*orsa::Vector(0,0,1);
        const double rotAxisOffset = acos(ref_rotAxis*rotAxis);
        
        std::vector< std::vector<double> > C, S, norm_C, norm_S;
        std::vector<double> J;
        orsa::convert(C, S, norm_C, norm_S, J,
                      paulMoment.get(), 
                      FromUnits(300,orsa::Unit::KM));
        std::vector<double> deltaMax;
        deltaMax.resize(degree+1); // [2] to [degree], but [0] and [1] are unused...
        double deltaRMS=0.0;
        unsigned int deltaRMS_count=0;
        for (unsigned int l=2; l<=degree; ++l) {
            deltaMax[l] = 0.0;
            for (unsigned int m=0; m<=l; ++m) {
                const double delta = fabs(norm_C[l][m]-ref_norm_C[l][m]);
                deltaRMS += orsa::square(delta);
                ++deltaRMS_count;
                if (delta>deltaMax[l]) {
                    deltaMax[l]=delta;
                }
                if (m!=0) {
                    const double delta = fabs(norm_S[l][m]-ref_norm_S[l][m]);
                    deltaRMS += orsa::square(delta);
                    ++deltaRMS_count;
                    if (delta>deltaMax[l]) {
                        deltaMax[l]=delta;
                    }
                }
            }
        }
        deltaRMS /= deltaRMS_count;
        deltaRMS = sqrt(deltaRMS);
        
        // output
        ORSA_DEBUG("SAMPLE: %8.1f %10.3e %10.3e %+8.1f %+8.1f %+8.1f %10.1f %10.1f %10.1f %8.1f %10.3e %10.3e %10.3e %10.3e %8.1f %10.6f",
                   (*val.coreDensity),
                   orsa::FromUnits(md->_coreVolume,orsa::Unit::KM,-3),
                   orsa::FromUnits(md->_coreVolume*val.coreDensity,orsa::Unit::KG,-1),
                   (*val.coreCenterX),
                   (*val.coreCenterY),
                   (*val.coreCenterZ),
                   (*val.coreRadiusX),
                   (*val.coreRadiusY),
                   (*val.coreRadiusZ),
                   md->_mantleDensity,
                   deltaMax[2],
                   deltaMax[3],
                   deltaMax[4],
                   deltaRMS,
                   centerOfMassOffset,
                   orsa::radToDeg()*rotAxisOffset);
        
    }    
    
    return 0;
}
