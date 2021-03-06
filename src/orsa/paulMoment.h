#ifndef _ORSA_PAUL_MOMENT_H_
#define _ORSA_PAUL_MOMENT_H_

#include <vector>
#include <string>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/double.h>
// #include <orsa/shape.h>
#include <orsa/vector.h>
// #include <orsa/massDistribution.h>
// #include <orsa/matrix.h>

namespace orsa {
  
    class PaulMoment : public osg::Referenced {
    public:
        PaulMoment(const unsigned int order);
    protected:
        virtual ~PaulMoment() { } 
    
    public:
        double M (const int i, const int j, const int k) const;
        double M_uncertainty (const int i, const int j, const int k) const;
    public:
        void setM (const double & val, const int i, const int j, const int k);
        void setM_uncertainty (const double & val, const int i, const int j, const int k);
    protected:
        std::vector< std::vector< std::vector<double> > > _M, _M_uncertainty;
    public:
        const unsigned int order;
    };
    
    bool translate(PaulMoment * const pm,
                   const PaulMoment * const pm0,
                   const orsa::Vector & delta);
    
    typedef std::vector< std::vector< std::vector<double> > >    triIndex_d;
    typedef std::vector< std::vector< std::vector<mpq_class> > > triIndex_mpq;
    const triIndex_mpq conversionCoefficients_C_integral(const size_t & l, const size_t & m);
    const triIndex_d   conversionCoefficients_C_plain(   const size_t & l, const size_t & m);
    const triIndex_d   conversionCoefficients_C_norm(    const size_t & l, const size_t & m);
    const triIndex_mpq conversionCoefficients_S_integral(const size_t & l, const size_t & m);
    const triIndex_d   conversionCoefficients_S_plain(   const size_t & l, const size_t & m);
    const triIndex_d   conversionCoefficients_S_norm(    const size_t & l, const size_t & m);
    
    // utility, just printing out values for now
    void convert(std::vector< std::vector<mpf_class> > & C,
                 std::vector< std::vector<mpf_class> > & S,
                 std::vector< std::vector<mpf_class> > & norm_C,
                 std::vector< std::vector<mpf_class> > & norm_S,
                 std::vector<mpf_class> & J,
                 const PaulMoment * const pm,
                 const double     & R0,
                 const bool         verbose=false);
    
    // tries to find a set of PaulMoments that matches the normalized C and S
    // the order of pm must be <= the order of norm_C and norm_S
    bool solve(PaulMoment * pm,
               const std::vector< std::vector<mpf_class> > & norm_C,
               const std::vector< std::vector<mpf_class> > & norm_S,
               const double     & R0);
    
    // a,b,c are the 3 semiaxes along x,y,z respectively
    void EllipsoidExpansion(PaulMoment   * pm,
                            const double & a,
                            const double & b,
                            const double & c);
    
}; // namespace orsa

#endif // _ORSA_PAUL_MOMENT_H_
