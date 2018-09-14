#ifndef _ORSA_CACHE_
#define _ORSA_CACHE_

#include <QQueue>

#define ORSA_CACHE_GET_CHECK 1

#if ORSA_CACHE_GET_CHECK
#include <orsa/debug.h>
#endif

#include <orsa/crash.h>

namespace orsa {
    
    //! Instead of Cache, this class should be called something like GetSetVariable...
    //
    template <typename T> class Cache {
    public:
        Cache() : _val(), _set(false), _locked(false) { }
    public:
        Cache(const T & val) : _val(val), _set(true), _locked(false) { }
    public:
        Cache(const Cache<T> & c) : _val(c._val), _set(c._set), _locked(c._locked) { }
    public:
        // virtual destructor is slow, and never needed so far...
        // virtual ~Cache() { }
        ~Cache() { }
        
#if ORSA_CACHE_GET_CHECK
    public:
        operator T const & () const {
            if (!_set) {
                ORSA_ERROR("returning unset value");
                orsa::crash();
            }
            return _val; 
        }
        const T & operator * () const {
            if (!_set) {
                ORSA_ERROR("returning unset value");
                orsa::crash();
            }
            return _val; 
        }       
        const T * operator -> () const {
            if (!_set) {
                ORSA_ERROR("returning unset value");
                orsa::crash();
            }
            return &_val; 
        }       
#else  
    public:
        operator T const & () const    { return  _val; }
        const T & operator * () const  { return  _val; }
        const T * operator -> () const { return &_val; }
#endif
        
    public:
        Cache<T> & operator = (const T & val) {
            if (canChange()) {
                _val = val;
                _set = true;
            }
            return (*this);
        }
    public:
        Cache<T> & operator = (const Cache<T> & c) {
            // ORSA_DEBUG("cache copy of cache");
            if (canChange()) {
                _val = c._val;
                _set = c._set;
                _locked = c._locked;
            }
            return (*this);
        }
    public:
        Cache<T> & operator += (const T & val) {
            if (canChange()) {
                _val += val;
            }
            return (*this);
        }
    public:
        Cache<T> & operator -= (const T & val) {
            if (canChange()) {
                _val -= val;
            }
            return (*this);
        }
    public:
        Cache<T> & operator *= (const T & val) {
            if (canChange()) {
                _val *= val;
            }
            return (*this);
        }
    public:
        Cache<T> & operator /= (const T & val) {
            if (canChange()) {
                _val /= val;
            }
            return (*this);
        }
        
    public:
        // set only if still not set
        bool conditionalSet(const T & val) {
            if (canChange()) {
                if (!_set) {
                    (*this) = val;
                    return true;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        }
        
    public:
        // set if not set, and
        // set if smaller than what is already set
        bool setIfSmaller(const T & val) {
            if (canChange()) {
                if (!_set) {
                    (*this) = val;
                    return true;
                } else {
                    if (val < _val) {
                        (*this) = val;
                        return true;
                    } else {
                        return false;
                    }
                }
            } else {
                return false;
            }
        }
        
    public:
        // set if not set, and
        // set if larger than what is already set
        bool setIfLarger(const T & val) {
            if (canChange()) {
                if (!_set) {
                    (*this) = val;
                    return true;
                } else {
                    if (val > _val) {
                        (*this) = val;
                        return true;
                    } else {
                        return false;
                    }
                }
            } else {
                return false;
            }
        }
        
    public:
        // use this to prevent changes to the variable
        // this has nothing to do with multi-threading locks...
        inline void lock()     const { _locked=true;   }
        inline void unlock()   const { _locked=false;  }
        inline bool isLocked() const { return _locked; }
    protected:
        inline bool canChange() const {
            if (_locked) {
                ORSA_ERROR("this variable is locked and cannot be changed");
                orsa::crash();
            }
            return (!_locked);
        }
        
    public:
        inline bool isSet() const { return _set; }
    public:
        void reset() { 
            if (_locked) {
                ORSA_ERROR("this variable is locked and cannot be changed");
                orsa::crash();
            }
            _set=false;
        }
    private:
        T    _val;
        bool _set;
        mutable bool _locked;
    };
    
    template <typename T> T operator + (const orsa::Cache<T> & c) {
        return T(c);
    }
    
    template <typename T> T operator - (const orsa::Cache<T> & c) {
        return -T(c);
    }
    
    template <typename T> T operator + (const orsa::Cache<T> & a, const orsa::Cache<T> & b) {
        return (T(a) + T(b));
    }
    
    template <typename T> T operator + (const orsa::Cache<T> & a, const T & b) {
        return (T(a) + b);
    }
    
    template <typename T> T operator + (const T & a, const orsa::Cache<T> & b) {
        return (a + T(b));
    }
    
    template <typename T> T operator - (const orsa::Cache<T> & a, const orsa::Cache<T> & b) {
        return (T(a) - T(b));
    }
    
    template <typename T> T operator - (const orsa::Cache<T> & a, const T & b) {
        return (T(a) - b);
    }
    
    template <typename T> T operator - (const T & a, const orsa::Cache<T> & b) {
        return (a - T(b));
    }
    
    template <typename T, typename P> P operator * (const orsa::Cache<T> & a, const orsa::Cache<T> & b) {
        return (T(a) * T(b));
    }
    
    template <typename T, typename P> P operator * (const orsa::Cache<T> & a, const T & b) {
        return (T(a) * b);
    }
    
    template <typename T, typename P> P operator * (const T & a, const orsa::Cache<T> & b) {
        return (a * T(b));
    }
    
    template <typename T, typename D> D operator / (const orsa::Cache<T> & a, const orsa::Cache<T> & b) {
        return (T(a) / T(b));
    }
    
    template <typename T, typename D> D operator / (const orsa::Cache<T> & a, const T & b) {
        return (T(a) / b);
    }
    
    template <typename T, typename D> D operator / (const T & a, const orsa::Cache<T> & b) {
        return (a / T(b));
    }
    
    template <typename T> bool operator == (const orsa::Cache<T> & a, const orsa::Cache<T> & b) {
        return (T(a) == T(b));
    }
    
    template <typename T> bool operator == (const orsa::Cache<T> & a, const T & b) {
        return (T(a) == b);
    }
    
    template <typename T> bool operator == (const T & a, const orsa::Cache<T> & b) {
        return (a == T(b));
    }
    
    template <typename T> bool operator != (const orsa::Cache<T> & a, const orsa::Cache<T> & b) {
        return (T(a) != T(b));
    }
    
    template <typename T> bool operator != (const orsa::Cache<T> & a, const T & b) {
        return (T(a) != b);
    }
    
    template <typename T> bool operator != (const T & a, const orsa::Cache<T> & b) {
        return (a != T(b));
    }
    
    template <typename T> bool operator <  (const orsa::Cache<T> & a, const orsa::Cache<T> & b) {
        return (T(a) <  T(b));
    }
    
    template <typename T> bool operator <  (const orsa::Cache<T> & a, const T & b) {
        return (T(a) <  b);
    }
    
    template <typename T> bool operator <  (const T & a, const orsa::Cache<T> & b) {
        return (a < T(b));
    }
    
    template <typename T> bool operator >  (const orsa::Cache<T> & a, const orsa::Cache<T> & b) {
        return (T(a) > T(b));
    }
    
    template <typename T> bool operator >  (const orsa::Cache<T> & a, const T & b) {
        return (T(a) >  b);
    }
    
    template <typename T> bool operator >  (const T & a, const orsa::Cache<T> & b) {
        return (a > T(b));
    }
    
    template <typename T> bool operator <= (const orsa::Cache<T> & a, const orsa::Cache<T> & b) {
        return (T(a) <= T(b));
    }
    
    template <typename T> bool operator <= (const orsa::Cache<T> & a, T & b) {
        return (T(a) <= b);
    }
    
    template <typename T> bool operator <= (const T & a, const orsa::Cache<T> & b) {
        return (a <= T(b));
    }
    
    template <typename T> bool operator >= (const orsa::Cache<T> & a, const orsa::Cache<T> & b) {
        return (T(a) >= T(b));
    }
    
    template <typename T> bool operator >= (const orsa::Cache<T> & a, const T & b) {
        return (T(a) >= b);
    }
    
    template <typename T> bool operator >= (const T & a, const orsa::Cache<T> & b) {
        return (a >= T(b));
    }
    
    
    
    template <class T> class CachedVars {
    public:
        CachedVars() {
            _s = 16;
        }
    public:
        CachedVars(const size_t max_size) {
            _s = max_size;
        }
    public:
        void setSize(const size_t max_size) {
            _s = max_size;
        }
    public:
        void insert(const T & val) const {
            _q.enqueue(val);
            while ((size_t)_q.size() > _s) {
                _q.dequeue();
            }
        }
    public:
        void clear() {
            _q.clear();
        }
    public:
        bool find(const T & val, T & found) const {
            qSort(_q.begin(),_q.end());
            typename QQueue<T>::const_iterator _it = qLowerBound(_q.begin(),_q.end(),val);
            if (_it == _q.end()) {
                ORSA_DEBUG("not found, size: %i _s: %i",_q.size(),_s);
                return false;
            }
            if ((*_it) == val) {
                ORSA_DEBUG(" -------------------- FOUND!!! ------------- size: %i _s: %i",_q.size(),_s);
                found = (*_it);
                return true;
            } else {
                ORSA_DEBUG("not found, size: %i _s: %i",_q.size(),_s);
                return false;
            }
        }
    protected:
        mutable QQueue<T> _q;
        size_t            _s;
    };
  
} // namespace orsa

#endif // _ORSA_CACHE_
