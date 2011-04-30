#include <orsaUtil/statisticalRanging.h>

#include <orsaUtil/multimin.h>

bool orsaUtil::statisticalRanging(orsaUtil::SR_AdaptiveIntervalVector & intervalVector,
                                  const double & initialThresholdLevel,
                                  const double & targetThresholdLevel,
                                  const double & minAdaptiveRange,
                                  const double & maxAdaptiveRange,
                                  const double & intervalResidualProbability,
                                  const size_t & targetSamples,
                                  const size_t & maxIter,
                                  // const orsaSolarSystem::OpticalObservationVector & allOpticalObs,
                                  const orsaUtil::SR_AuxiliaryData * aux) {
    
    std::vector< orsa::Cache<double> > vecMinRange;
    vecMinRange.resize(aux->vecSize);
    std::vector< orsa::Cache<double> > vecMaxRange;
    vecMaxRange.resize(aux->vecSize);
    {
        // test: determine the maximum range for each observation
        osg::ref_ptr<orsaUtil::MultiminMinMaxRange> mmr = new orsaUtil::MultiminMinMaxRange;
        
        // run on ALL j and k because range search is asymmetric (determines at t=tj)
        for (unsigned int j=0; j<aux->vecSize; ++j) {
            for (unsigned int k=0; k<aux->vecSize; ++k) {
                if (j==k) continue;
                orsa::Cache<double> minRange, maxRange;
                mmr->getMinMaxRange(minRange, maxRange,
                                    minAdaptiveRange,
                                    maxAdaptiveRange,
                                    aux->R_s[j],aux->R_o[j],aux->u_o2a[j],aux->vecRMS[j],aux->allOpticalObs[j]->epoch,
                                    aux->R_s[k],aux->R_o[k],aux->u_o2a[k],aux->vecRMS[k],aux->allOpticalObs[k]->epoch,
                                    orsaSolarSystem::Data::GMSun(),
                                    3.0);
#warning nominal max error here... 3.0? 2.5?
                
                // maxRange is only for first set of arguments == [j]
                // vecMaxRange[j].setIfLarger(maxRange);
                // if (maxRange > 0.0) {
                if (minRange.isSet()) {
                    if (minRange > 0.0) {
                        vecMinRange[j].setIfSmaller(minRange);
                    }
                }
                if (maxRange.isSet()) {
                    if (maxRange > 0.0) {
                        vecMaxRange[j].setIfLarger(maxRange);
                    }
                }
                
                /* if (maxRange.isSet()) {
                   ORSA_DEBUG("max range distance for obs [%02i] and [%02i] = %8.3g [AU]",
                   j,
                   k,
                   orsa::FromUnits(maxRange,orsa::Unit::AU,-1));
                   }
                */
                
            }
            // ORSA_DEBUG("set: %i",vecMaxRange[j].isSet());
            if (vecMinRange[j].isSet()) {
                ORSA_DEBUG("min range distance for obs [%02i] = %8.3f [AU]",
                           j,
                           orsa::FromUnits(vecMinRange[j],orsa::Unit::AU,-1));
            } else {
#warning is it correct to force it here?
                vecMinRange[j] = minAdaptiveRange; 
            }
            if (vecMaxRange[j].isSet()) {
                ORSA_DEBUG("max range distance for obs [%02i] = %8.3f [AU]",
                           j,
                           orsa::FromUnits(vecMaxRange[j],orsa::Unit::AU,-1));
            } else {
#warning is it correct to force it here?
                vecMaxRange[j] = maxAdaptiveRange; 
            }
        }           
    }
    
    intervalVector.resize(aux->vecSize);
    for (size_t j=0; j<aux->vecSize; ++j) {
        // intervalVector[j] = new SR_AdaptiveInterval(0.0, 1.0, 1.0e-6, 100.0, 1.0, 1000);
        intervalVector[j] = new SR_AdaptiveInterval(std::max(minAdaptiveRange,(*vecMinRange[j])),
                                                    std::min((*vecMaxRange[j]),maxAdaptiveRange),
                                                    intervalResidualProbability, // "1-confidence level" for this interval, different from the chisq-level
                                                    initialThresholdLevel,
                                                    targetThresholdLevel,
                                                    targetSamples,
                                                    aux);
    }
    
    osg::ref_ptr<SR_AdaptiveMonteCarlo> mc = new SR_AdaptiveMonteCarlo;
    
    mc->run(intervalVector,maxIter);
    
    return true;
}
