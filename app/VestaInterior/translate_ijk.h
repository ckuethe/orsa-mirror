#ifndef _TRANSLATE_IJK_
#define _TRANSLATE_IJK_

void translate(std::vector< std::vector< std::vector<double> > > & tN,
               const std::vector< std::vector< std::vector<double> > > & N,
               const orsa::Vector & delta) {
    const size_t degree = N.size()-1;
    tN.resize(degree+1);
    for (size_t ni=0; ni<=degree; ++ni) {
        tN[ni].resize(degree+1-ni);
        for (size_t nj=0; nj<=degree-ni; ++nj) {
            tN[ni][nj].resize(degree+1-ni-nj);
            for (size_t nk=0; nk<=degree-ni-nj; ++nk) {
                tN[ni][nj][nk] = 0.0;
            }
        }
    }
    for (size_t ni=0; ni<=degree; ++ni) {
        for (size_t nj=0; nj<=degree-ni; ++nj) {
            for (size_t nk=0; nk<=degree-ni-nj; ++nk) {
                for (size_t bi=0; bi<=ni; ++bi) {
                    for (size_t bj=0; bj<=nj; ++bj) {
                        for (size_t bk=0; bk<=nk; ++bk) {
                            tN[ni][nj][nk] +=
                                mpz_class(orsa::binomial(ni,bi) *
                                          orsa::binomial(nj,bj) *
                                          orsa::binomial(nk,bk)).get_d() *
                                orsa::int_pow(delta.getX(),bi) *
                                orsa::int_pow(delta.getY(),bj) *
                                orsa::int_pow(delta.getZ(),bk) *
                                N[ni-bi][nj-bj][nk-bk];                            
                        }
                    }
                }
            }
        }
    }
}

#endif // _TRANSLATE_IJK_
