#ifndef _MOVING_LAYERS_MULTIFIT_H_
#define _MOVING_LAYERS_MULTIFIT_H_

#include <orsa/multifit.h>

#include "CubicChebyshevMassDistribution.h"

typedef dd_real simplex_T;

class MovingLayersMultifit : public orsa::Multifit {
public:
    MovingLayersMultifit() :
        orsa::Multifit()
        { }
public:
    class SHLayerData {
    public:
        size_t degree;
        double volume;
        double excessDensity;
        orsa::Vector v0;
    };
public:
    // all input vars
    size_t SH_degree;
    size_t SH_size;
    LayerData::EllipsoidLayerVectorType ellipsoidLayerVector; // this is never changed
    std::vector<SHLayerData> shLayerData;
    std::vector< std::vector<mpf_class> > uniformShape_norm_C;
    std::vector< std::vector<mpf_class> > uniformShape_norm_S;
    mpf_class                             uniformShape_IzzMR2;
    double R0_plate;
    double R0_gravity;
    double bulkDensity;
    double layersTotalMassFraction;
    double uniformShapeMassFraction;
    orsa::Vector sampled_CM;
    osg::ref_ptr<SimplexIntegration<simplex_T> > si;
    osg::ref_ptr<orsaPDS::RadioScienceGravityData> gravityData;
    
protected:
    void computeAllFunctionCalls(const orsa::MultifitParameters *, 
                                 const orsa::MultifitData       *,
                                 const computeAllCallsMode) const;
public:
    double fun(const orsa::MultifitParameters * par, 
               const orsa::MultifitData       * data,
               const unsigned int p, 
               const int          d,
               const unsigned int row) const;
protected:
    void save_CCMDF(const orsa::MultifitParameters *) const;
    void singleIterationDone(const orsa::MultifitParameters *) const;
    void success(const orsa::MultifitParameters *) const;
    void runCompleted(const bool /* success */, const orsa::MultifitParameters *) const;
protected:
    int f_gsl(const gsl_vector * parameters, 
              void             * dataPoints, 
              gsl_vector       * f) {
        
        int retval =  orsa::Multifit::f_gsl(parameters, 
                                            dataPoints, 
                                            f);
        
        const orsa::MultifitData * data = (orsa::MultifitData *) dataPoints;
        for (unsigned int j=0; j<data->size(); ++j) {
            ORSA_DEBUG("f[%02i] = %12.6g", j, gsl_vector_get(f,j));
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
                ORSA_DEBUG("df[%02i][%02i] = %12.6g", j, k, gsl_matrix_get(J,j,k));
            }
        }
        
        return retval;
    }
  
};

#endif // _MOVING_LAYERS_MULTIFIT_H_
