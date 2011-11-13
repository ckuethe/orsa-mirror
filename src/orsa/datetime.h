#ifndef _ORSA_DATETIME_H_
#define _ORSA_DATETIME_H_

#include <orsa/double.h>
#include <orsa/debug.h>
#include <orsa/unit.h>
#include <string>

namespace orsa {

    class Time {
    public:
        Time() { }
    public:
        Time(const Time & t) : _mu_sec(t._mu_sec) { }
    public:
        Time(const mpz_class mu_sec) : _mu_sec(mu_sec) { }
    public:
        Time(const mpz_class d,
             const mpz_class H,
             const mpz_class M,
             const mpz_class S,
             const mpz_class mu_S) {
            _mu_sec = mu_S + 1000000 *
                (S + 60 *
                 (M + 60 *
                  (H + 24 * d)));
        }
    public:
        /* inline Time & operator = (const mpz_class & rhs) {
           if (rhs != 0) {
           ORSA_DEBUG("warning: setting time from non-zero integer = %Zi",rhs.get_mpz_t());
           }
           _mu_sec = rhs;
           return * this;
           }
        */
    public:
        ~Time() { }
    public:
        bool set(const mpz_class d,
                 const mpz_class H,
                 const mpz_class M,
                 const mpz_class S,
                 const mpz_class mu_S) {
            _mu_sec = mu_S + 1000000 *
                (S + 60 *
                 (M + 60 *
                  (H + 24 * d)));
            return true;
        }
    public:
        const mpz_class & getMuSec() const {
            return _mu_sec;
        }
    public:
        double get_d() const {
            return (FromUnits(_mu_sec.get_d(),Unit::MICROSECOND));
        }
    public:
        inline Time & operator += (const Time & rhs) {
            _mu_sec += rhs._mu_sec;
            return * this;
        }
    public:
        inline Time & operator -= (const Time & rhs) {
            _mu_sec -= rhs._mu_sec;
            return * this;
        }
    public:
        inline Time & operator *= (const mpz_class & rhs) {
            _mu_sec *= rhs;
            return * this;
        }
    
    public:
        inline const Time operator + (const Time & rhs) const {
            Time _t(*this);
            _t += rhs;
            return _t;
        }
    public:
        inline const Time operator - (const Time & rhs) const {
            Time _t(*this);
            _t -= rhs;
            return _t;
        }
    
    public:
        inline const Time operator / (const mpz_class & rhs) const {
            Time _t(*this);
            _t._mu_sec /= rhs;
            return _t;
        }

    public:
        inline const Time operator + () const {
            Time _t(*this);
            return _t;
        }
    public:
        inline const Time operator - () const {
            Time _t(0);
            _t -= (*this);
            return _t;
        }

    public:
        inline bool operator < (const Time & rhs) const {
            return (_mu_sec < rhs._mu_sec);
        }
    public:
        inline bool operator <= (const Time & rhs) const {
            return (_mu_sec <= rhs._mu_sec);
        }
        inline bool operator > (const Time & rhs) const {
            return (_mu_sec > rhs._mu_sec);
        }
    public:
        inline bool operator >= (const Time & rhs) const {
            return (_mu_sec >= rhs._mu_sec);
        }
    public:
        inline bool operator == (const Time & rhs) const {
            return (_mu_sec == rhs._mu_sec);
        }
    public:
        inline bool operator != (const Time & rhs) const {
            return (_mu_sec != rhs._mu_sec);
        }
    
    protected:
        mpz_class _mu_sec;
    
    // Added for long long int
    protected:
        mpz_class _tmp;
    
    public:
        Time(long long int d,
             long long int H,
             long long int M,
             long long int S,
             long long int mu_S) {
             
             double test = ((double)(mu_S))+((double)(1000000.))*(((double)(S))+((double)(60.))*(((double)(M))+((double)(60.))*(((double)(H))+((double)(24.))*((double)(d)))));
            
            if (fabs(test)<((double)(9223372036854775800LL))) {
                long long int musec = ((long long int)(mu_S))+((long long int)(1000000))*(((long long int)(S))+((long long int)(60))*(((long long int)(M))+((long long int)(60))*(((long long int)(H))+((long long int)(24))*((long long int)(d)))));
                _mu_sec.set_str(lltostr(musec), 10);
            } else {
                mpz_class mpzd;
                mpzd.set_str(lltostr(d), 10);

                mpz_class mpzH;
                mpzH.set_str(lltostr(H), 10);

                mpz_class mpzM;
                mpzM.set_str(lltostr(M), 10);

                mpz_class mpzS;
                mpzS.set_str(lltostr(S), 10);

                mpz_class mpzmu_S;
                mpzmu_S.set_str(lltostr(mu_S), 10);
                
                _mu_sec = mpzmu_S + 1000000 *
                    (mpzS + 60 *
                     (mpzM + 60 *
                      (mpzH + 24 * mpzd)));
            }
        }
        
    public:
        bool set(long long int d,
                 long long int H,
                 long long int M,
                 long long int S,
                 long long int mu_S) {
             
             double test = ((double)(mu_S))+((double)(1000000.))*(((double)(S))+((double)(60.))*(((double)(M))+((double)(60.))*(((double)(H))+((double)(24.))*((double)(d)))));
            
            if (fabs(test)<((double)(9223372036854775800LL))) {
                long long int musec = ((long long int)(mu_S))+((long long int)(1000000))*(((long long int)(S))+((long long int)(60))*(((long long int)(M))+((long long int)(60))*(((long long int)(H))+((long long int)(24))*((long long int)(d)))));
                _mu_sec.set_str(lltostr(musec), 10);
            } else {
                mpz_class mpzd;
                mpzd.set_str(lltostr(d), 10);

                mpz_class mpzH;
                mpzH.set_str(lltostr(H), 10);

                mpz_class mpzM;
                mpzM.set_str(lltostr(M), 10);

                mpz_class mpzS;
                mpzS.set_str(lltostr(S), 10);

                mpz_class mpzmu_S;
                mpzmu_S.set_str(lltostr(mu_S), 10);
                
                _mu_sec = mpzmu_S + 1000000 *
                    (mpzS + 60 *
                     (mpzM + 60 *
                      (mpzH + 24 * mpzd)));
            }
            return true;
        }

    public:
        inline Time& operator=(const Time & rhs)
        {
            _mu_sec = rhs._mu_sec;
            return *this;
        }
    
    public:
        long long int get_ll(bool &ok) const {
            long long int result = 0;
            ok = false;
            double dble = _mu_sec.get_d();
            if (fabs(dble)<((double)(9223372036854775800LL))) { // 2^63-1 = 9223372036854775800
                result = std::stoll(_mu_sec.get_str(), NULL, 10);
                ok = true;
            }
            return result;
        }

    public:
        inline Time & operator += (long long int musec) {
            _tmp.set_str(lltostr(musec), 10);
            _mu_sec += _tmp;
            return * this;
        }
        
    public:
        inline Time & operator -= (long long int musec) {
            _tmp.set_str(lltostr(musec), 10);
            _mu_sec -= _tmp;
            return * this;
        }
    public:
        inline const Time operator + (long long int musec) const {
            Time _t(*this);
            _t += musec;
            return _t;
        }
    public:
        inline const Time operator - (long long int musec) const {
            Time _t(*this);
            _t -= musec;
            return _t;
        }
        
    public:
        inline bool operator < (long long int musec) {
            _tmp.set_str(lltostr(musec), 10);
            return (_mu_sec < _tmp);
        }
    public:
        inline bool operator <= (long long int musec) {
            _tmp.set_str(lltostr(musec), 10);
            return (_mu_sec <= _tmp);
        }
        inline bool operator > (long long int musec) {
            _tmp.set_str(lltostr(musec), 10);
            return (_mu_sec > _tmp);
        }
    public:
        inline bool operator >= (long long int musec) {
            _tmp.set_str(lltostr(musec), 10);
            return (_mu_sec >= _tmp);
        }
    public:
        inline bool operator == (long long int musec) {
            _tmp.set_str(lltostr(musec), 10);
            return (_mu_sec == _tmp);
        }
    public:
        inline bool operator != (long long int musec) {
            _tmp.set_str(lltostr(musec), 10);
            return (_mu_sec != _tmp);
        }
    
    public:
        static inline void strreverse(char* begin, char* end)
        {
            char aux;
            while (end > begin)
                aux = *end, *end-- = *begin, *begin++ = aux;
        }
    public:
        static inline std::string lltostr(long long int value)
        {
            return std::to_string(value);
        }
    };

    inline const Time operator * (const Time      & lhs,
                                  const mpz_class & rhs) {
        return Time(lhs.getMuSec()*rhs);
    }

    inline const Time operator * (const mpz_class & lhs,
                                  const Time      & rhs) {
        return Time(lhs*rhs.getMuSec());
    }

} // namespace orsa

#endif // _ORSA_DATETIME_H_
