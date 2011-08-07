#ifndef _ORSA_DOUBLE_
#define _ORSA_DOUBLE_

#include <cstdio>
#include <gmpxx.h>

#include <complex>

#include <orsa/debug.h>
#include <iostream>
#include <string>
#include <cstdlib>

namespace orsa {
    
    // computes i! (factorial)
    mpz_class factorial(const mpz_class & i,
                        const bool & cache=true);
    
    // computes i!! (bi-factorial)
    mpz_class bi_factorial(const mpz_class & i,
                           const bool & cache=true);
    
    mpz_class binomial(const mpz_class & n,
                       const mpz_class & k,
                       const bool & cache=true);
    
    // (-1)^l
    int power_sign(const mpz_class & l);
    
    template <class T> T int_pow(const T & x, const int & p);
    
    template <class T> T square(const T & x) { return (x*x); }
    
    template <class T> T cube(const T & x) { return (x*x*x); }
    
    const double & epsilon();
    
    const double & pi();
    
    const double & halfpi();
    
    const double & twopi();
    
    const double & pisquared();
    
    const double & radToDeg();
    
    const double & degToRad();
    
    const double & radToArcmin();
    
    const double & arcminToRad();
    
    const double & radToArcsec();
    
    const double & arcsecToRad();
    
    int kronecker(const mpz_class & i,
                  const mpz_class & j);
    
    mpf_class pochhammer(const mpf_class & x, 
                         const mpz_class & n);
    
    mpz_class pochhammer(const mpz_class & a, 
                         const mpz_class & n);
    
    // call orsa::crash() if x is not finite (+-Inf or NaN), useful for debugging
    void check(const double & x);
    
} // namespace orsa

#endif // _ORSA_DOUBLE_
