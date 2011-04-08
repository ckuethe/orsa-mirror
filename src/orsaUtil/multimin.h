#ifndef _ORSA_UTIL_MULTIMIN_
#define _ORSA_UTIL_MULTIMIN_

#include <orsa/datetime.h>
#include <orsa/multimin.h>
#include <orsa/vector.h>

namespace orsaUtil {
    
    // multimin orbital velocity
    
    class MultiminOrbitalVelocity : public orsa::Multimin {
    public:
        double fun(const orsa::MultiminParameters * par) const;
    public:
        orsa::Vector getOrbitalVelocity(
            const orsa::Vector & R1_,
            const orsa::Time   & t1_,
            const orsa::Vector & R2_,
            const orsa::Time   & t2_,
            const double       mu_);
    protected:
        // utility
        inline orsa::Vector getVel(const orsa::MultiminParameters * par) const {
            return (par->get("Vr")*u_R+par->get("Vn")*u_N);
        }
    protected:
        orsa::Cache<orsa::Vector> R1, R2;
        orsa::Cache<double> mu;
        orsa::Cache<double> dt;
        orsa::Cache<orsa::Vector> u_R, u_N, u_L; // R=Radial; N=Normal to {L,R}; L=Angular momentum
    };

    // multimin maximum range

    class MultiminMaximumRange : public orsa::Multimin {
    public:
        double fun(const orsa::MultiminParameters * par) const;
    public:
        double getMaximumRange(const orsa::Vector & R_s_1_,
                               const orsa::Vector & R_o_1_,
                               const orsa::Vector & u_o2a_1_,
                               const double       & sigma_1_arcsec_,
                               const orsa::Time   & epoch_1_,
                               const orsa::Vector & R_s_2_,
                               const orsa::Vector & R_o_2_,
                               const orsa::Vector & u_o2a_2_,
                               const double       & sigma_2_arcsec_,
                               const orsa::Time   & epoch_2_,
                               const double       & GM_s_,
                               const double       & sigma_factor_);
    protected:
        orsa::Cache<orsa::Vector> R_s_1, R_o_1, u_o2a_1;
        orsa::Cache<double>       sigma_1_arcsec;
        orsa::Cache<orsa::Time>   epoch_1;
    protected:
        orsa::Cache<orsa::Vector> R_s_2, R_o_2, u_o2a_2;
        orsa::Cache<double>       sigma_2_arcsec;
        orsa::Cache<orsa::Time>   epoch_2;
    protected:
        orsa::Cache<double>       GM_s;
        orsa::Cache<double>       sigma_factor;
    };
    
}; // namespace orsaUtil

#endif // _ORSA_UTIL_MULTIMIN_
