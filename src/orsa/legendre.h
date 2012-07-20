#ifndef _ORSA_LEGENDRE_
#define _ORSA_LEGENDRE_

#include <orsa/cache.h>
#include <orsa/double.h>

#include <vector>

namespace orsa {
    
    inline mpf_class LegendreP(const mpz_class & l,
                               const mpz_class & m,
                               const mpf_class & x) {
        mpf_class sum_plus=0;
        mpf_class sum_minus=0;
        // for(mpz_class k=0; k<=l/2; ++k) {
        for (mpz_class k=l/2; k>=0; --k) {
            /* ORSA_DEBUG("k: %Zi   sum_plus: %Fg   sum_minus: %Fg",
               k.get_mpz_t(),
               sum_plus.get_mpf_t(),
               sum_minus.get_mpf_t());
            */
            
            if (k%2==0) sum_plus +=
                // power_sign(k) *
                binomial(l,k) *
                binomial(2*l-2*k,l) *
                pochhammer(mpz_class(l-m-2*k+1),m) *
                int_pow(x,l-m-2*k);
            else sum_minus +=
                // power_sign(k) *
                binomial(l,k) *
                binomial(2*l-2*k,l) *
                pochhammer(mpz_class(l-m-2*k+1),m) *
                int_pow(x,l-m-2*k);
        }
        mpf_class sum = sum_plus - sum_minus;
        // ORSA_DEBUG("sum: %Fg",sum.get_mpf_t());
        
        sum *= pow(mpf_class(1-x*x).get_d(),m.get_d()/2.0);
        // ORSA_DEBUG("sum: %Fg",sum.get_mpf_t());
        sum /= int_pow(mpz_class(2),l);
        // ORSA_DEBUG("sum: %Fg",sum.get_mpf_t());
        
        return sum;
    }
    
} // namespace orsa

#endif // _ORSA_LEGENDRE_
