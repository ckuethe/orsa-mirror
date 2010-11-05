#ifndef _ORSA_BINOMIAL_H_
#define _ORSA_BINOMIAL_H_

double singleProbabilityTerm(const double & p,
                             const double & N,
                             const double & R,
                             const bool & cache=true) {
    return (p * orsa::binomial(N,R,cache).get_d() * orsa::int_pow(p,R) * orsa::int_pow(1-p,N-R));
}       

double cumulativeProbability(const double & p,
                             const double & N,
                             const double & R,
                             const bool & cache=true) {
    double runningProb = 0.0;
    for (unsigned int Z=R; Z<=N; ++Z) {
        runningProb += singleProbabilityTerm(p,Z,R,cache);
    }
    return runningProb;
}

bool bisectionSearch(unsigned int & N_minus,
                     unsigned int & N_plus,
                     const unsigned int & R,
                     const double & p,
                     const double & probThreshold,
                     const bool & cache=true) {
    if (N_minus > N_plus) return false;
    double p_minus =  cumulativeProbability(p, N_minus, R, cache);
    double p_plus  =  cumulativeProbability(p, N_plus,  R, cache);
    if (p_minus > probThreshold) {
        if (N_minus==R) {
            N_plus = R+1;
            return true;
        } else {
            ORSA_DEBUG("problem: N_minus is above threshold");
            return false;
        }
    }
    if (p_plus  < probThreshold) {
        ORSA_DEBUG("problem: N_plus is below threshold");
        return false;
    }
    unsigned int N_middle;
    double p_middle;
    do {
        N_middle = (N_plus+N_minus)/2;
        if (N_middle==N_minus) break;
        if (N_middle==N_plus)  break;
        p_middle = cumulativeProbability(p, N_middle, R, cache);
        ORSA_DEBUG("N: %i  cP: %g",N_middle,p_middle);
        if (p_middle <= probThreshold) {
            N_minus = N_middle;
            p_minus = p_middle;
        } else {
            N_plus = N_middle;
            p_plus = p_middle;
        } 
    } while (N_plus-N_minus>1);
    return true;
}

bool BinomialConfidenceInterval(unsigned int & N_low,
                                unsigned int & N_mean,
                                unsigned int & N_high,
                                const double & CL /* confidence level (between 0.0 and 1.0) */,
                                const unsigned int & R /* found objects */,
                                const double & p /* probability */,
                                const bool & cache=true) {
    // check input
    if ((p<=0.0) || (p>1.0)) {
        ORSA_DEBUG("invalid probability value");
        return false;
    } else if ((CL<=0.0) || (CL>=1.0)) {
        ORSA_DEBUG("invalid confidence level value");
        return false;
    }
    
    N_mean = lrint(((R+1)/p)-1);
    
    // notice how the R==0 case is handled in a special way for both N_low and N_high
    
    if (R==0) {
        N_low = R;
    } else {
        // search N_low
        const double probThreshold = (1.0-CL)/2;
        unsigned int N_minus = R;
        unsigned int N_plus  = N_mean;
        if (bisectionSearch(N_minus,N_plus,R,p,probThreshold,cache)) {
            N_low = N_minus;
        } else {
            return false;
        }
    }
    
    {
        // search N_high
        const double probThreshold = (R==0 ? CL : (1.0+CL)/2);
        unsigned int N_minus = N_mean;
        unsigned int N_plus  = 5*N_mean;
        if (bisectionSearch(N_minus,N_plus,R,p,probThreshold,cache)) {
            N_high = N_plus;
        } else {
            return false;
        }
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
