#ifndef _ORSA_BINOMIAL_H_
#define _ORSA_BINOMIAL_H_

double singleProbabilityTerm(const double & p,
                             const unsigned int & N,
                             const unsigned int & R,
                             const bool & cache=true) {
    double sPT = (p * orsa::binomial(N,R,cache).get_d() * orsa::int_pow(p,R) * orsa::int_pow(1-p,N-R));
    if (!finite(sPT)) {
        // need MPF
        sPT = mpf_class(mpf_class(p) * orsa::binomial(N,R,cache) * orsa::int_pow(mpf_class(p),R) * orsa::int_pow(mpf_class(1-p),N-R)).get_d();
    }
    if (!finite(sPT)) {
        ORSA_DEBUG("PROBLEMS -- p: %g  N: %i  R: %i   binomial(N,R): %Zi (float=%g)",
                   p,N,R,
                   orsa::binomial(N,R,cache).get_mpz_t(),
                   orsa::binomial(N,R,cache).get_d());
        orsa::crash();
    }
    return sPT;
}       

double cumulativeProbability(const double & p,
                             const unsigned int & N,
                             const unsigned int & R,
                             const bool & cache=true) {
    double runningProb = 0.0;
    for (unsigned int Z=R; Z<=N; ++Z) {
        const double sPT = singleProbabilityTerm(p,Z,R,cache);
        /* if (!finite(sPT)) {
           ORSA_DEBUG("Z=%i --> adding sPT=%g to rP=%g",Z,sPT,runningProb);
           }
        */
        runningProb += sPT;
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
            N_plus = R;
            return true;
        } else {
            ORSA_DEBUG("problem: N_minus=%i is above threshold, retrying...",N_minus);
            --N_minus;
            return bisectionSearch(N_minus,
                                   N_plus,
                                   R,
                                   p,
                                   probThreshold,
                                   cache);
        }
    }
    if (p_plus < probThreshold) {
        ORSA_DEBUG("problem: N_plus=%i is below threshold",N_plus);
        return false;
    }
    unsigned int N_middle;
    double p_middle;
    do {
        N_middle = (N_plus+N_minus)/2;
        if (N_middle==N_minus) break;
        if (N_middle==N_plus)  break;
        p_middle = cumulativeProbability(p, N_middle, R, cache);
        // ORSA_DEBUG("N: %i  cP: %g",N_middle,p_middle);
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
                                const double & CL,      /* confidence level (between 0.0 and 1.0) */
                                const unsigned int & R, /* found objects */
                                const double & p,       /* cumulative probability */
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
        // ORSA_DEBUG("searching N_low...");
        const double probThreshold = (1.0-CL)/2;
        unsigned int N_minus = R;
        unsigned int N_plus  = N_mean+1;
        if (bisectionSearch(N_minus,N_plus,R,p,probThreshold,cache)) {
            N_low = N_minus;
        } else {
            return false;
        }
    }
    
    {
        // search N_high
        // ORSA_DEBUG("searching N_high...");
        const double probThreshold = (R==0 ? CL : (1.0+CL)/2);
        unsigned int N_minus = N_mean;
        unsigned int N_plus  = 5*(N_mean+1);
        if (bisectionSearch(N_minus,N_plus,R,p,probThreshold,cache)) {
            N_high = N_plus;
        } else {
            return false;
        }
    }
    
    /* ORSA_DEBUG("R: %i   p: %g   N_low: %i   N_mean = %g -> %i   N_high: %i   (C.L.: %g\%)",
       R,
       p,
       N_low,
       ((R+1)/p)-1,
       N_mean,
       N_high,
       100*CL);
    */
    
    return true;
}

#endif // _ORSA_BINOMIAL_H_
