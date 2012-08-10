#ifndef _MOVING_LAYERS_MULTIFIT_H_
#define _MOVING_LAYERS_MULTIFIT_H_

#include <orsa/multifit.h>

class MovingLayersMultifit : public orsa::Multifit {
public:
    MovingLayersMultifit() :
        orsa::Multifit()
        { }
protected:
    void singleIterationDone(const orsa::MultifitParameters *) const {
        ORSA_DEBUG("--MARK--");
    }
public:
    double fun(const orsa::MultifitParameters * par, 
               const orsa::MultifitData       * data,
               const unsigned int p, 
               const int          d,
               const unsigned int row) const {
        
    }
protected:
    int f_gsl(const gsl_vector * parameters, 
              void             * dataPoints, 
              gsl_vector       * f) {
        
        int retval =  orsa::Multifit::f_gsl(parameters, 
                                            dataPoints, 
                                            f);
        
        const orsa::MultifitData * data = (orsa::MultifitData *) dataPoints;
        for (unsigned int j=0; j<data->size(); ++j) {
            ORSA_DEBUG("f[%02i] = %10.3f", j, gsl_vector_get(f,j));
        }
        
        return retval;
    }
protected:
    int df_gsl (const gsl_vector * v, 
                void             * dataPoints, 
                gsl_matrix       * J) {
    
        int retval = Multifit::df_gsl(v, 
                                      dataPoints, 
                                      J);
    
        const orsa::MultifitData * data = (orsa::MultifitData *) dataPoints;
        for (unsigned int k=0; k<_par->totalSize(); ++k) {
            for (unsigned int j=0; j<data->size(); ++j) {
                ORSA_DEBUG("df[%02i][%02i] = %10.3f", j, k, gsl_matrix_get(J,j,k));
            }
        }
        
        return retval;
    }
  
};

#endif // _MOVING_LAYERS_MULTIFIT_H_
