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
  // parent_ has to be set before calling CheckInput_ otherwise one gets a
  // segfault!
  parent_ = parent;
  value_ = value;

  CheckInput_(value);
  
  time_.reserve(value.size());
  for (int i = 0; i < value.size(); i++) {
    time_[i] = i;
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
FlexibleInput<T>::FlexibleInput(cyclus::Agent* parent,
                                std::vector<T> value, 
                                std::vector<int> time) {
  // parent_ has to be set before calling CheckInput_ otherwise one gets a
  // segfault!
  value_ = value;
  time_ = time;
  parent_ = parent;

  CheckInput_(value, time);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
T FlexibleInput<T>::UpdateValue() {
  // Get current time with t = 0 being the entrance of parent_ in the simulation
  int t = parent_->context()->time() - parent_->enter_time();
  
  if (t >= time_[current_idx_] && t < time_[current_idx_+1]) {
    return value_[current_idx_];
  } else if (t == time_[current_idx_+1]) {
    ++current_idx_;
    return value_[current_idx_];
  } else {
    std::stringstream ss;
    ss << "Agent '" << parent_->prototype()  << "' of spec '" << parent_->spec() 
       << "' with enter_time '" << parent_->enter_time() 
       << "' has passed the invalid timestamp '" << t << "' at time '" 
       << parent_->context()->time() << "' to a FlexibleInput variable.\n";
    
    throw cyclus::ValueError(ss.str());
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
void FlexibleInput<T>::CheckInput_(const std::vector<T>& value) {
  int lifetime = parent_->lifetime();
  if (lifetime == -1) {
    lifetime = parent_->context()->sim_info().duration;
  }

  if (value.size() > lifetime) {
    std::stringstream ss;
    ss << "While initialising agent '" << parent_->prototype()  
       << "' of spec '" << parent_->spec() << "' at time'" 
       << parent_->context()->time() << "' a problem appeared:\n"
       << "The value vector passed to FlexibleInput contains too many "
       << "elements.\n";
    
    throw cyclus::ValueError(ss.str());
 }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
void FlexibleInput<T>::CheckInput_(const std::vector<T>& value,
                                   const std::vector<int>& time_vec) {
  CheckInput_(value);
  if (value.size() != time_vec.size()) {
    std::stringstream ss;
    ss << "While initialising agent '" << parent_->prototype()  
       << "' of spec '" << parent_->spec() << "' at time '" 
       << parent_->context()->time() << "' a problem appeared:\n"
       << "time and value vectors do not have the same size.\n"
       << "size of time vector: " << time_vec.size()
       << ", size of value vector: " << value.size() << "\n";
    
    throw cyclus::ValueError(ss.str());
  }

}

}  // namespace misoenrichment
