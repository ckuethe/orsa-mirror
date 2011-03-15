#include <orsaSPICE/spiceBodyTranslationalCallback.h>

#include <orsaSPICE/spice.h>

using namespace orsa;
using namespace orsaSPICE;

SpiceBodyTranslationalCallback::SpiceBodyTranslationalCallback(const std::string & name) : 
    orsa::PrecomputedTranslationalBodyProperty(),
    _name(name) { }


SpiceBodyTranslationalCallback::SpiceBodyTranslationalCallback(const SpiceBodyTranslationalCallback & sbtc) :
    orsa::PrecomputedTranslationalBodyProperty(),
    _name(sbtc._name) { 
    if (sbtc._previousTime.isSet()) {
        _position     = sbtc._position;
        _velocity     = sbtc._velocity;
        _previousTime = sbtc._previousTime;
    }
}

orsa::Vector SpiceBodyTranslationalCallback::position() const { return _position; }

orsa::Vector SpiceBodyTranslationalCallback::velocity() const { return _velocity; }

bool SpiceBodyTranslationalCallback::update(const orsa::Time & t) {
  
    if (_previousTime.isSet()) {
        if (_previousTime == t) {
            // ORSA_DEBUG("cached...");
            return true;
        }    
    }
  
    _previousTime = t;
  
    orsa::Vector r, v;
  
    SPICE::instance()->getPosVel(_name,
                                 t,
                                 r,
                                 v);
    _position = r;
    _velocity = v;
  
    return true;
}

void  SpiceBodyTranslationalCallback::lock() {
    mutex.lock();
}

void SpiceBodyTranslationalCallback::unlock() {
    mutex.unlock();
}
