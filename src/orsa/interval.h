#ifndef _ORSA_INTERVAL_
#define _ORSA_INTERVAL_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <algorithm>
#include <deque>

// #include <QList>

#include <orsa/cache.h>
#include <orsa/crash.h>
#include <orsa/debug.h>

namespace orsa {
    
    // NOTE: it should not be possible to modify an Interval from outside....
    // for the moment, if you modify it, you MUST call update() when you're done! 
    
    // stored data: only unique values...
    template <typename T> class Interval : public osg::Referenced {
    public:
        typedef std::deque<T> DataType;
    public:
        Interval() : Referenced(true) { 
            _store_data = false;
        }
    protected:
        virtual ~Interval() { }
    public:
        bool insert(const T & val, const bool onlyIfExtending, const bool replaceIfDouble) {
            if (_store_data) {
                if (size() == 0) {
                    _data.push_front(val);
                } else {
                    if (val < _min) {
                        _data.push_front(val);
                    } else if (val > _max) {
                        _data.push_back(val);
                    } else if (!onlyIfExtending) {
                        //
                        typename DataType::iterator _it = lower_bound(_data.begin(),_data.end(),val);
                        if ((*_it) != val) {
                            _data.insert(_it,val);
                        } else {
                            if (replaceIfDouble) {
                                *_it = val;
                            } else {
                                ORSA_DEBUG("called Interval<T>::insert() with duplicate entry");
                                return false;
                            }
                        }
                    } else {
                        // ORSA_DEBUG("insert failed, onlyIfExtending=%i",onlyIfExtending);
                        return false;
                    }
                }
                update();
                // ORSA_DEBUG("size: %i",size());
            } else {
                if (_min.isSet()) {
                    if (val < _min) {
                        _min = val;
                    }
                } else {
                    _min = val;
                }
                if (_max.isSet()) {
                    if (val > _max) {
                        _max = val;
                    }
                } else {
                    _max = val;
                }
            }
            return true;
        }
    public:
        bool remove(const T & val) {
            if (_store_data) {
                typename DataType::iterator _it = lower_bound(_data.begin(),_data.end(),val);
                if ((*_it) == val) {
                    _it = _data.erase(_it);
                    update();
                }
                return true;
            } else {
                return false;
            }
        }
        
    public:
        bool reset() { 
            _data.clear();
            _min.reset();
            _max.reset();
            // ORSA_DEBUG("called reset(), this = %x",this);
            return true;
        }
    public:
        bool update() {
            if (_store_data && size()) {
                _min = *(_data.begin());
                _max = *(--_data.end());
            }  
            return true;
        }
    
    public:
        bool enableDataStoring() {
            _store_data = true;
            return true;
        }
    public:
        bool disableDataStoring() {
            _store_data = false;
            _data.clear();
            return true;
        }
    public:    
        bool isStoringData() const {
            return _store_data;
        }
    
    public:
        const T & min() const {
            return _min;
        }
    
    public:
        const T & max() const {
            return _max;
        }
    
    public:
        bool getSubInterval(const T & val, T & sub_min, T & sub_max) const {
            
            // ORSA_DEBUG("size: %i   this = %x",size(),this);
            
            if (size()==0) {
                ORSA_DEBUG("returning false...   size: %i   this = %x",size(),this);
                // orsa::crash();
                return false;
            }
            
            if (!isStoringData()) {
                ORSA_ERROR("this interval is not storing data");
                return false;
            }
            // ORSA_DEBUG("//1//");
            if (val == max()) {
                sub_min = sub_max = max();
                // ORSA_DEBUG("max...");
                return true;
            }
            // ORSA_DEBUG("//2//");
            if (val == min()) {
                sub_min = sub_max = min();
                // ORSA_DEBUG("min...");
                return true;
            }
            // ORSA_DEBUG("//3//");
            if (val < min()) {
                ORSA_ERROR("value requested is below minimum, size: %i",size());
                sub_min = sub_max = min();
                return false;
            }
            // ORSA_DEBUG("//4//");
            if (val > max()) {
                ORSA_ERROR("value requested is above maximum, size: %i",size());
                sub_min = sub_max = max();
                return false;
            }
            // ORSA_DEBUG("//5//");
            
            if (size() == 1) {
                sub_max = max();
                sub_min = min();
                return true;
            }
            
            if (size() == 2) {
                sub_max = max();
                sub_min = min();
                return true;
            }
            
            // ORSA_DEBUG("//6//");
            //
            // ORSA_DEBUG("getSubInterval(...) is searching...");
            // std::cerr << "[s]";
            //
            typename DataType::const_iterator _it = lower_bound(_data.begin(),_data.end(),val);
            // typename DataType::const_iterator _it = qLowerBound(_data.begin(),_data.end(),val);
            // typename DataType::const_iterator _it = _num_edge_const_lower_bound(val);
            //
            if ((*_it) == val) {
	
                sub_max = sub_min = (*_it);
                
                // ORSA_DEBUG("//7//");
                
                return true;
                /* 
                // IF you need some code in here, the interval has some problems (i.e. was modified from outside the class, without calling update() when done)
                } else if ((*_it) > max()) {
                ORSA_DEBUG("//XXX// (got it!) [1]");
                sub_min = sub_max = max();
                return true;
                } else if ((*_it) == min()) {
                ORSA_DEBUG("//XXX// (got it!) [2]");
                sub_min = sub_max = min();
                return true;
                */
            } else {
                sub_max = (*_it);
                sub_min = (*(--_it));
                
                // ORSA_DEBUG("//8//");
                
                return true;	  
            }	  
            
            ORSA_DEBUG("returning false...   size: %i   this = %x",size(),this);
            
            return false;
        }	
        
    public:
        size_t size() const {
            return _data.size();
        }
    
    protected:
        bool _store_data;
    public:
        const DataType & getData() const {
            return _data;
        }
    public:
        // non-const version too
        // you MUST call update() after any change...
        DataType & getData() {
            return _data;
        }
    protected:
        DataType _data;
    protected:
        Cache<T> _min, _max;
    };
  
}; // namespace orsa

#endif // _ORSA_INTERVAL_
