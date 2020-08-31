#include "flexible_input.h"

#include <algorithm>
#include <sstream>

#include "context.h"
#include "error.h"

namespace misoenrichment {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
FlexibleInput<T>::FlexibleInput() {;}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
FlexibleInput<T>::FlexibleInput(cyclus::Agent* parent, 
                                std::vector<T> value) {
  value_ = value;
  time_.reserve(value.size());
  for (int i = 0; i < value.size(); i++) {
    time_[i] = i;
  }
  time_it_ = time_.begin();
  CheckInput_(parent, value);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
FlexibleInput<T>::FlexibleInput(cyclus::Agent* parent,
                                std::vector<T> value, 
                                std::vector<int> time) {
  value_ = value;
  time_ = time;
  time_it_ = time_.begin();
  CheckInput_(parent, value, time);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
T FlexibleInput<T>::UpdateValue(cyclus::Agent* parent) {
  // Get current time with t = 0 being the entrance of parent in the 
  // simulation.
  int t = parent->context()->time() - parent->enter_time();

  // The second conditional takes the ending of the time vector into
  // account. If the last element is reached, *(time_it_+1) is not 
  // evaluated.
  if (t >= *time_it_ && (time_it_+1 == time_.end() || t < *(time_it_+1))) {
    return value_[time_it_ - time_.begin()];
  } else if (t == *(time_it_+1)) {
    ++time_it_;
    return value_[time_it_ - time_.begin()];
  } else {
    std::stringstream ss;
    ss << "Agent '" << parent->prototype()  
       << "' of spec '" << parent->spec() 
       << "' with enter_time '" << parent->enter_time()
       << "' has passed the invalid timestamp '" << t << "' at time '" 
       << parent->context()->time() << "' to a FlexibleInput variable.\n";
    
    throw cyclus::ValueError(ss.str());
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
void FlexibleInput<T>::CheckInput_(cyclus::Agent* parent, 
                                   const std::vector<T>& value) {
  int lifetime = parent->lifetime();
  if (lifetime == -1) {
    lifetime = parent->context()->sim_info().duration;
  }

  if (value.size() > lifetime) {
    std::stringstream ss;
    ss << "While initialising agent '" << parent->prototype()  
       << "' of spec '" << parent->spec() << "' at time "
       << parent->context()->time() << " a problem appeared:\n"
       << "The value vector passed to FlexibleInput contains too many "
       << "elements.\n";
    
   throw cyclus::ValueError(ss.str());
 }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
void FlexibleInput<T>::CheckInput_(cyclus::Agent* parent, 
                                   const std::vector<T>& value,
                                   const std::vector<int>& time_vec) {
  CheckInput_(parent, value);
  if (value.size() != time_vec.size()) {
    std::stringstream ss;
    ss << "While initialising agent '" << parent->prototype()  
       << "' of spec '" << parent->spec() << "' at time '" 
       << parent->context()->time() << "' a problem appeared:\n"
       << "time and value vectors do not have the same size.\n"
       << "size of time vector: " << time_vec.size()
       << ", size of value vector: " << value.size() << "\n";
    
    throw cyclus::ValueError(ss.str());
  }

  if (time_vec[0] != 0) {
    std::stringstream ss;
    ss << "While initialising agent '" << parent->prototype()  
       << "' of spec '" << parent->spec() << "' at time '" 
       << parent->context()->time() << "' a problem appeared:\n"
       << "the first element of the time vector must be '0' (initial "
       << "value).\n";
    
    throw cyclus::ValueError(ss.str());
  }
}

}  // namespace misoenrichment
