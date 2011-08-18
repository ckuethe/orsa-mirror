#ifndef _ORSA_CHEBYSHEV_
#define _ORSA_CHEBYSHEV_

#include <vector>

#include <orsa/double.h>
#include <orsa/debug.h>

namespace orsa {
    
    inline void ChebyshevT(std::vector<double> & T,
                           const size_t & n,
                           const double & x) {
        // ORSA_DEBUG("n: %2i x: %+6.3f",n,x);
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
    
    inline const std::vector<mpz_class> & ChebyshevTcoeff(const size_t & n) {
        static std::vector< std::vector<mpz_class> > coeff;
        if (coeff.size() > n) {
            return coeff[n];
        } else {
            const size_t oldSize = coeff.size();
            coeff.resize(n+1);
            for (size_t k=oldSize; k<=n; ++k) {
                coeff[k].resize(k+1);
                if (k==0) {
                    coeff[0][0] = 1;
                } else if (k==1) {
                    coeff[1][0] = 0;
                    coeff[1][1] = 1;
                } else {
                    for (size_t j=0; j<=k; ++j) {
                        coeff[k][j] = 0;
                        if (j>0) coeff[k][j] += 2*coeff[k-1][j-1];
                        if (j<=k-2) coeff[k][j] -= coeff[k-2][j];
                    }
                }
            }
            return coeff[n];
        }
    }
    
    inline const std::vector<double> & ChebyshevTextrema(const size_t & n) {
        static std::vector< std::vector<double> > ext;
        if (ext.size() > n) {
            return ext[n];
        } else {
            const size_t oldSize = ext.size();
            ext.resize(n+1);
            for (size_t k=oldSize; k<=n; ++k) {
                ext[k].resize(k+1);
                if (k==0) {
                    ext[0][0] = 1.0; // T_0=1 flat, so just pick a point
                } else {
                    for (size_t j=0; j<=k; ++j) {
                        ext[k][j] = cos(orsa::pi()*j/k);
                    }
                }
            }
            return ext[n];
        }
    }
    
} // namespace orsa

#endif // _ORSA_CHEBYSHEV_
