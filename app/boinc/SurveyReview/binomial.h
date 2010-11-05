#ifndef _ORSA_BINOMIAL_H_
#define _ORSA_BINOMIAL_H_

bool BinomialConfidenceInterval(unsigned int & N_low,
                                unsigned int & N_mean,
                                unsigned int & N_high,
                                const double & CL /* confidence level (between 0.0 and 1.0) */,
                                const unsigned int & R /* found objects */,
                                const double & p /* probability */) {
    // check input
    if ((p<=0.0) || (p>1.0)) {
        ORSA_DEBUG("invalid probability value");
        return false;
    } else if ((CL<=0.0) || (CL>1.0)) {
        ORSA_DEBUG("invalid confidence level value");
        return false;
    }
    
    N_mean = lrint(((R+1)/p)-1);
    
    {
        // search N_low
        const double probThreshold = (1.0-CL)/2;
        double runningBinProb = 0.0;
        unsigned int Z=R;
        while (1) {
            runningBinProb +=
                p *
                orsa::binomial(Z,R,false).get_d() *
                orsa::int_pow(p,R) *
                orsa::int_pow(1-p,Z-R);
            ORSA_DEBUG("(low) N: %i  rBP: %g",Z,runningBinProb);
            if (runningBinProb >= probThreshold) {
                if (Z > R) {
                    --Z;
                }
                break;
            } else {
                ++Z;
            }
        }
        N_low = Z;
    }
    
    {
        // search N_high
        const double probThreshold = (1.0+CL)/2;
        double runningBinProb = 0.0;
        unsigned int Z=R;
        while (1) {
            runningBinProb +=
                p *
                orsa::binomial(Z,R,false).get_d() *
                orsa::int_pow(p,R) *
                orsa::int_pow(1-p,Z-R);
            ORSA_DEBUG("(high) N: %i  rBP: %g",Z,runningBinProb);
            if (runningBinProb >= probThreshold) {
                break;
            } else {
                ++Z;
            }
        }
        N_high = Z;
    }
    
    ORSA_DEBUG("R: %i   p: %g   N_low: %i   N_mean = %g -> %i   N_high: %i   (C.L.: %g\%)",
               R,
               p,
               N_low,
               ((R+1)/p)-1,
               N_mean,
               N_high,
               100*CL);
    
    return true;
}

#endif // _ORSA_BINOMIAL_H_
