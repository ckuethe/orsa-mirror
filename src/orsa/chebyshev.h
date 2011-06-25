#ifndef _ORSA_CHEBYSHEV_
#define _ORSA_CHEBYSHEV_

#include <vector>

namespace orsa {
    
    inline void ChebyshevT(std::vector<double> & T,
                           const size_t & n,
                           const double & x) {
        T.resize(n+1);
        T[0] = 1;
        if (n==0) return;
        T[1] = x;
        if (n==1) return;
        for (size_t s=2; s<=n; ++s) {
            T[s] = 2*x*T[s-1]-T[s-2];
        }
    }
    
    inline double ChebyshevT(const size_t & n, const double & x) {
        std::vector<double> T;
        orsa::ChebyshevT(T,n,x);
        const double retVal = T[n];
        return retVal;
    }
    
} // namespace orsa

#endif _ORSA_CHEBYSHEV_
