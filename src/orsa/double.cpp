#include <orsa/double.h>
#include <orsa/cache.h>

#include <vector>

using namespace orsa;

mpz_class orsa::factorial(const mpz_class & i,
                          const bool & cache) {
    // ORSA_DEBUG("i: %Zi  cache: %i",i.get_mpz_t(),cache);
    if (cache) {
        static std::vector< orsa::Cache<mpz_class> > _factorial_table;
        const unsigned long int index = i.get_ui();
        // static mpz_class _mpz_one("1");
        if (i <= 1) {
            return 1;
        } else if (_factorial_table.size() > i) {
            if (!_factorial_table[index].isSet()) {
                _factorial_table[index] = i*factorial(i-1,cache);
            }
            return _factorial_table[index];
        } else {
            _factorial_table.resize(index+1);
            _factorial_table[index] = i*factorial(i-1,cache);
            return _factorial_table[index];
        }
    } else {
        if (i <= 1) {
            return 1;
        }
        return (i*factorial(i-1,cache));    
    }   
}

mpz_class orsa::bi_factorial(const mpz_class & i,
                             const bool & cache) {
    if (cache) {
        // ORSA_DEBUG("bf: %Zi",i.get_mpz_t());
        static std::vector< orsa::Cache<mpz_class> > _bi_factorial_table;
        const unsigned long int index = i.get_ui();
        // static mpz_class _mpz_one("1");
        // static mpz_class _mpz_two("2");
        if (i <= 1) {
            return 1;
        } else if (_bi_factorial_table.size() > i) {
            if (!_bi_factorial_table[index].isSet()) {
                _bi_factorial_table[index] = i*bi_factorial(i-2,cache);
            }
            return _bi_factorial_table[index];
        } else {
            _bi_factorial_table.resize(index+1);
            _bi_factorial_table[index] = i*bi_factorial(i-2,cache);
            return _bi_factorial_table[index];
        }
    } else {
        if (i <= mpz_class("1")) {
            return mpz_class("1");
        }
        return (i*bi_factorial(i-2,cache));
    }
}

mpz_class orsa::binomial(const mpz_class & n,
                         const mpz_class & k,
                         const bool & cache) {
    // ORSA_DEBUG("n: %Zi   k: %Zi   cache: %i",n.get_mpz_t(),k.get_mpz_t(),cache);
    if ((k<0) || (k>n)) return 0;
    return factorial(n,cache)/(factorial(k,cache)*factorial(n-k,cache));
}

int orsa::power_sign(const mpz_class & l) {
    if ((((l%2)+2)%2)==1) {
        return -1;
    } else {
        return  1;
    }
}

/* double orsa::int_pow(const double & x, 
   const int    & p) {
   if (p ==  2) return x*x;
   if (p ==  1) return x;
   if (p ==  0) return 1;
   if (p == -1) return 1/x;
   double _pow = x;
   const int max_k = abs(p);
   for (int k=1; k<max_k; ++k) {
   _pow *= x;
   }
   if (p < 0) _pow = 1/_pow;
   return _pow;
   }
*/

template <class T> T orsa::int_pow(const T & x,
                                   const int & p) {
    if (p ==  2) return x*x;
    if (p ==  1) return x;
    if (p ==  0) return 1;
    if (p == -1) return 1/x;
    T _pow = x;
    const int max_k = abs(p);
    for (int k=1; k<max_k; ++k) {
        _pow *= x;
    }
    if (p < 0) _pow = T(1)/_pow;
    return _pow;
}

template <class T> T orsa::int_pow(const T & x,
                                   const mpz_class & p) {
    if (p ==  2) return x*x;
    if (p ==  1) return x;
    if (p ==  0) return 1;
    if (p == -1) return 1/x;
    T _pow = x;
    const mpz_class max_k = abs(p);
    for (mpz_class k=1; k<max_k; ++k) {
        _pow *= x;
    }
    if (p < 0) _pow = T(1)/_pow;
    return _pow;
}

template int orsa::int_pow(const int & x,
                           const int & p);

template double orsa::int_pow(const double & x,
                              const int & p);

template mpf_class orsa::int_pow(const mpf_class & x,
                                 const int & p);

template mpz_class orsa::int_pow(const mpz_class & z,
                                 const int & p);

template int orsa::int_pow(const int & x,
                           const mpz_class & p);

template mpz_class orsa::int_pow(const mpz_class & z,
                                 const mpz_class & p);

template mpf_class orsa::int_pow(const mpf_class & x,
                                 const mpz_class & p);

const double & orsa::epsilon() {
    static const double _eps = __DBL_EPSILON__; /* 2.2204460492503131e-16 */
    return _eps;
}

const double & orsa::pi() {
    static const double _pi = 3.14159265358979323846;
    return _pi;
}

const double & orsa::halfpi() {
    static const double _halfpi = 0.5*orsa::pi();
    return _halfpi;
}

const double & orsa::twopi() {
    static const double _twopi = 2*orsa::pi();
    return _twopi;
}

const double & orsa::pisquared() {
    static const double _pisquared = orsa::pi()*orsa::pi();
    return _pisquared;
}

const double & orsa::radToDeg() {
    static const double _radToDeg = 180/orsa::pi();
    return _radToDeg;
}

const double & orsa::degToRad() {
    static const double _degToRad = orsa::pi()/180;
    return _degToRad;
}

const double & orsa::radToArcmin() {
    static const double _radToArcmin = 180*60/orsa::pi();
    return _radToArcmin;
}

const double & orsa::arcminToRad() {
    static const double _arcminToRad = orsa::pi()/(180*60);
    return _arcminToRad;
}

const double & orsa::radToArcsec() {
    static const double _radToArcsec = 180*3600/orsa::pi();
    return _radToArcsec;
}

const double & orsa::arcsecToRad() {
    static const double _arcsecToRad = orsa::pi()/(180*3600);
    return _arcsecToRad;
}

int orsa::kronecker(const mpz_class & i,
                    const mpz_class & j) {
    if (i==j) {
        return 1;
    } else {
        return 0;
    }
}

mpf_class orsa::pochhammer(const mpf_class & x,
                           const mpz_class & n) {
    if (n == 0) {
        return 1;
    } else if (n == 1) {
        return x;
    } else if (n > 0) {
        mpf_class _result("1");
        mpz_class _k("0");
        while (_k<n) {
            _result *= (x+_k);
            ++_k;
        }
        return _result;
    } else {
        mpf_class _result("1");
        mpz_class _k("-1");
        while (_k>=n) {
            // ORSA_DEBUG("k: %Zi n: %Zi",_k.get_mpz_t(),n.get_mpz_t());
            _result /= (x+_k);
            --_k;
        }      
        return _result;
    }
}

mpz_class orsa::pochhammer(const mpz_class & a,
                           const mpz_class & n) {
    if (n < 0) {
        ORSA_DEBUG("must use the mpf version of pochhammer for negative n");
        return 0;
    } else if (n == 0) {
        return 1;
    } else if (n == 1) {
        return a;
    } else {
        mpz_class _result("1");
        mpz_class _k("0");
        while (_k < n) { 
            // ORSA_DEBUG("k: %Zi n: %Zi",_k.get_mpz_t(),n.get_mpz_t());
            _result *= (a+_k);
            ++_k;
        }
        return _result;
    }
}

void orsa::check(const double & x) {
    if (!finite(x)) {
        ORSA_DEBUG("detected a double with %g value",x);
        orsa::crash();
    }
}
