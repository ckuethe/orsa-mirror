#include <orsa/box.h>

#include <osg/BoundingBox>

using namespace orsa;

Box::Box() { }

Box::Box(const double & xMin,
         const double & xMax,
         const double & yMin,
         const double & yMax,
         const double & zMin,
         const double & zMax) {
    set(xMin, xMax, yMin, yMax, zMin, zMax);
}

Box::Box(const orsa::Vector & v1,
         const orsa::Vector & v2) {
    set(v1, v2);
}

void Box::set(const double & xMin,
              const double & xMax,
              const double & yMin,
              const double & yMax,
              const double & zMin,
              const double & zMax) {
    _xMin = xMin;
    _xMax = xMax;
    _yMin = yMin;
    _yMax = yMax;
    _zMin = zMin;
    _zMax = zMax;
}

void Box::set(const orsa::Vector & v1,
              const orsa::Vector & v2) {
    if (v1.getX() < v2.getX()) {
        _xMin = v1.getX();
        _xMax = v2.getX();
    } else {
        _xMin = v2.getX();
        _xMax = v1.getX();
    }
  
    if (v1.getY() < v2.getY()) {
        _yMin = v1.getY();
        _yMax = v2.getY();
    } else {
        _yMin = v2.getY();
        _yMax = v1.getY();
    }
  
    if (v1.getZ() < v2.getZ()) {
        _zMin = v1.getZ();
        _zMax = v2.getZ();
    } else {
        _zMin = v2.getZ();
        _zMax = v1.getZ();
    }
}

bool Box::isSet() const {
    return (_xMin.isSet() &&
            _xMax.isSet() &&
            _yMin.isSet() &&
            _yMax.isSet() &&
            _zMin.isSet() &&
            _zMax.isSet());
}

void Box::reset() {
    _xMin.reset();
    _xMax.reset();
    _yMin.reset();
    _yMax.reset();
    _zMin.reset();
    _zMax.reset();
}

double Box::volume() const {
    return fabs((_xMax-_xMin)*
                (_yMax-_yMin)*
                (_zMax-_zMin));
}

bool Box::isInside(const orsa::Vector & v) const {
    if (v.getX() < _xMin) return false;
    if (v.getX() > _xMax) return false;
    if (v.getY() < _yMin) return false;
    if (v.getY() > _yMax) return false;
    if (v.getZ() < _zMin) return false;
    if (v.getZ() > _zMax) return false;
    return true;
}

osg::BoundingBox Box::getOSGBoundingBox() const {
    return osg::BoundingBox(_xMin,
                            _yMin,
                            _zMin,
                            _xMax,
                            _yMax,
                            _zMax);
}

