#ifndef _SHAPE_H_
#define _SHAPE_H_

// a single point to define the type of shape model input file type

#include <orsa/shape.h>

#include "vesta.h"
#include "gaskell.h"
#include "eros_shape.h"

// pick one here
class InputShape : public ErosShape { }; // PDS, Bennu,...
// class InputShape : public VestaShape { };
// class InputShape : public GaskellPlateModel { };

#endif // _SHAPE_H_
