#ifndef _ROTATE_IJK_
#define _ROTATE_IJK_

#include <orsa/double.h>

void rotate(std::vector< std::vector< std::vector<double> > > & rN,
            const std::vector< std::vector< std::vector<double> > > & N,
            const orsa::Matrix & m) {
    const size_t degree = N.size()-1;
    rN.resize(degree+1);
    for (size_t ni=0; ni<=degree; ++ni) {
        rN[ni].resize(degree+1-ni);
        for (size_t nj=0; nj<=degree-ni; ++nj) {
            rN[ni][nj].resize(degree+1-ni-nj);
            for (size_t nk=0; nk<=degree-ni-nj; ++nk) {
                rN[ni][nj][nk] = 0.0;
            }
        }
    }
    const double & m11 = m.getM11();
    const double & m12 = m.getM12();
    const double & m13 = m.getM13();
    const double & m21 = m.getM21();
    const double & m22 = m.getM22();
    const double & m23 = m.getM23();
    const double & m31 = m.getM31();
    const double & m32 = m.getM32();
    const double & m33 = m.getM33();
    for (size_t ni=0; ni<=degree; ++ni) {
        for (size_t nj=0; nj<=degree-ni; ++nj) {
            for (size_t nk=0; nk<=degree-ni-nj; ++nk) {
                for (size_t s11=0; s11<=ni; ++s11) {
                    for (size_t s12=0; s12<=ni; ++s12) {
                        for (size_t s13=0; s13<=ni; ++s13) {
                            if (s11+s12+s13==ni) {
                                for (size_t s21=0; s21<=nj; ++s21) {
                                    for (size_t s22=0; s22<=nj; ++s22) {
                                        for (size_t s23=0; s23<=nj; ++s23) {
                                            if (s21+s22+s23==nj) {
                                                for (size_t s31=0; s31<=nk; ++s31) {
                                                    for (size_t s32=0; s32<=nk; ++s32) {
                                                        for (size_t s33=0; s33<=nk; ++s33) {
                                                            if (s31+s32+s33==nk) {
                                                                rN[ni][nj][nk] +=
                                                                    mpz_class(orsa::factorial(ni)/(orsa::factorial(s11)*orsa::factorial(s12)*orsa::factorial(s13)) *
                                                                              orsa::factorial(nj)/(orsa::factorial(s21)*orsa::factorial(s22)*orsa::factorial(s23)) *
                                                                              orsa::factorial(nk)/(orsa::factorial(s31)*orsa::factorial(s32)*orsa::factorial(s33))).get_d() *
                                                                    orsa::int_pow(m11,s11) *
                                                                    orsa::int_pow(m12,s12) *
                                                                    orsa::int_pow(m13,s13) *
                                                                    orsa::int_pow(m21,s21) *
                                                                    orsa::int_pow(m22,s22) *
                                                                    orsa::int_pow(m23,s23) *
                                                                    orsa::int_pow(m31,s31) *
                                                                    orsa::int_pow(m32,s32) *
                                                                    orsa::int_pow(m33,s33) *
                                                                    N[s11+s21+s31][s12+s22+s32][s13+s23+s33];            
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#endif // _ROTATE_IJK_
