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
        Cache() : _val(), _set(false) { }
    public:
        Cache(const T & val) : _val(val), _set(true) { }
    public:
        Cache(const Cache<T> & c) : _val(c._val), _set(c._set) { }
    public:
        // virtual destructor is slow, and never needed so far...
        // virtual ~Cache() { }
        ~Cache() { }
    public:
#if ORSA_CACHE_GET_CHECK
    public:
        /* operator T() const {
           if (!_set) {
           ORSA_ERROR("returning unset value");
           orsa::crash();
           }
           return _val; 
           }
        */
        operator T const & () const {
            if (!_set) {
                ORSA_ERROR("returning unset value");
                orsa::crash();
            }
            return _val; 
        }
        /* T get() const { 
           if (!_set) {
           ORSA_ERROR("returning unset value");
           orsa::crash();
           }
           return _val; 
           }
           const T & getRef() const {
           if (!_set) {
           ORSA_ERROR("returning unset value");
           orsa::crash();
           }
           return _val;
           } 
           const T * getPtr() const {
           if (!_set) {
           ORSA_ERROR("returning unset value");
           orsa::crash();
           }
           return & _val;
           }
        */
#else  
    public:
        // inline operator T() const { return _val; }
        inline operator T const & () const { return _val; }
        /* inline T get() const { return _val; }
           inline const T & getRef() const { return _val; }
           inline const T * getPtr() const { return & _val; }
        */
#endif
    protected:
        /* bool set(const T & val) { 
           _val=val; 
           _set=true; 
           return _set;
           }
        */
    public:
        Cache<T> & operator = (const T & val) {
            // set(val);
            _val = val;
            _set = true;
            return (*this);
        }
    public:
        Cache<T> & operator += (const T & val) {
            _val += val;
            return (*this);
        }
    public:
        Cache<T> & operator -= (const T & val) {
            _val -= val;
            return (*this);
        }


        // some useful methods... could also be implemented using T(Cache<T>).fun()
    public:
        // useful when T is a std::string
        const char * c_str ( ) const { return _val.c_str(); }
    public:
        // useful when T is a orsa::Time
        double get_d() const { return _val.get_d(); }
    public:
        // useful when T is a orsa::Angle
        const double & getRad() const { return _val.getRad(); }
        
        
    public:
        // set only if still not set
        bool conditionalSet(const T & val) {
            if (!_set) {
                (*this) = val;
                return true;
            } else {
                return false;
            }
        }
    public:
        inline bool isSet() const { return _set; }
    public:
        void reset() { 
            _set=false;
        }
    private:
        T    _val;
        bool _set;
    };
    
    template <typename T> T operator + (const orsa::Cache<T> & c) {
        T t = c;
        return t;
    }
    
    template <typename T> T operator - (const orsa::Cache<T> & c) {
        T t = c;
        return -t;
    }

    template <typename T> T operator + (const orsa::Cache<T> & a, const orsa::Cache<T> & b) {
        T t = a;
        t += b;
        return t;
    }
    
    template <typename T> T operator + (const orsa::Cache<T> & a, const T & b) {
        T t = a;
        t += b;
        return t;
    }
    
    template <typename T> T operator + (const T & a, const orsa::Cache<T> & b) {
        T t = a;
        t += b;
        return t;
    }
    
    template <typename T> T operator - (const orsa::Cache<T> & a, const orsa::Cache<T> & b) {
        T t = a;
        t -= b;
        return t;
    }
    
    template <typename T> T operator - (const orsa::Cache<T> & a, const T & b) {
        T t = a;
        t -= b;
        return t;
    }
    
    template <typename T> T operator - (const T & a, const orsa::Cache<T> & b) {
        T t = a;
        t -= b;
        return t;
    }
    
    template <typename T, typename P> P operator * (const orsa::Cache<T> & a, const orsa::Cache<T> & b) {
        return (T(a) *  T(b));
    }
    
    template <typename T, typename P> P operator * (const orsa::Cache<T> & a, const T & b) {
        return (T(a) *  b);
    }
    
    template <typename T, typename P> P operator * (const T & a, const orsa::Cache<T> & b) {
        return (a *  T(b));
    }
    
    template <typename T, typename D> D operator / (const orsa::Cache<T> & a, const orsa::Cache<T> & b) {
        return (T(a) /  T(b));
    }
    
    template <typename T, typename D> D operator / (const orsa::Cache<T> & a, const T & b) {
        return (T(a) /  b);
    }
    
    template <typename T, typename D> D operator / (const T & a, const orsa::Cache<T> & b) {
        return (a /  T(b));
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
        return (a <  T(b));
    }
    
    template <typename T> bool operator >  (const orsa::Cache<T> & a, const orsa::Cache<T> & b) {
        return (T(a) >  T(b));
    }
    
    template <typename T> bool operator >  (const orsa::Cache<T> & a, const T & b) {
        return (T(a) >  b);
    }
    
    template <typename T> bool operator >  (const T & a, const orsa::Cache<T> & b) {
        return (a >  T(b));
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
