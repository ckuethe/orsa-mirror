#include <orsa/util.h>
#include <orsaInputOutput/MPC_asteroid.h>
#include <orsaInputOutput/MPC_observations.h>
#include <orsaSolarSystem/datetime.h>

int main(int argc, char ** argv) {
    
    // outputs the discovery observations for a given obsCode, along with the orbital info of the object
    
    orsa::Debug::instance()->initTimer();
  
    if (argc != 3) {
        printf("Usage: %s <obsFile> <obsCode>\n",argv[0]);
        exit(0);
    }
    
    // ORSA_DEBUG("process ID: %i",getpid());
    
    osg::ref_ptr<orsaInputOutput::MPCObservationsFile> obsFile = new orsaInputOutput::MPCObservationsFile;
    obsFile->setFileName(argv[1]);
    obsFile->select_obsCode = argv[2];
    obsFile->select_discovery = true;
    obsFile->read();
    
    osg::ref_ptr<orsaInputOutput::MPCAsteroidFile> orbitFile = new orsaInputOutput::MPCAsteroidFile;
    orbitFile->setFileName("NEA.txt");
    orbitFile->read();
    
    // match each observation with an object and its orbit
    for (unsigned int korb=0; korb<orbitFile->_data.size(); ++korb) {
        orsaInputOutput::MPCAsteroidDataElement asteroid;
        osg::ref_ptr<orsaSolarSystem::OpticalObservation> observation;
        bool found=false;
        for (unsigned int kobs=0; kobs<obsFile->_data.size(); ++kobs) {
            osg::ref_ptr<orsaSolarSystem::OpticalObservation> obs = dynamic_cast<orsaSolarSystem::OpticalObservation *> (obsFile->_data[kobs].get());
            if (obs) {
                if (obs->designation.isSet() && orbitFile->_data[korb].designation.isSet()) {
                    if (obs->designation.getRef() == orbitFile->_data[korb].designation.getRef()) {
                        asteroid    = orbitFile->_data[korb];
                        observation = obs;
                        found       = true;
                        break;
                    }
                }
                if (obs->number.isSet() && orbitFile->_data[korb].number.isSet()) {
                    if (obs->number.getRef() == orbitFile->_data[korb].number.getRef()) {
                        asteroid    = orbitFile->_data[korb];
                        observation = obs;
                        found       = true;
                        break;
                    }
                }
            }
        }
        if (found) {
            ORSA_DEBUG("%f %12.6f %5.2f %10.6f %10.6f %10.6f %9i %12s %s",
                       orsaSolarSystem::timeToJulian(observation->epoch.getRef()),
                       orsaSolarSystem::fractionalYear(observation->epoch.getRef()),
                       asteroid.H.getRef(),
                       orsa::FromUnits(asteroid.orbit.getRef().a,orsa::Unit::AU,-1),
                       asteroid.orbit.getRef().e,
                       asteroid.orbit.getRef().i*orsa::radToDeg(),
                       (asteroid.number.isSet()?asteroid.number.getRef():0),
                       (asteroid.designation.isSet()?asteroid.designation.getRef().c_str():"nodes"),
                       observation->obsCode.getRef().c_str());
        }
    }
    
    return 0;
}
