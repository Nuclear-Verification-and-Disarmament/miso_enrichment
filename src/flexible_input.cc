#include "flexible_input.h"

#include <algorithm>
#include <sstream>

#include "error.h"

namespace misoenrichment {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
FlexibleInput<T>::FlexibleInput(cyclus::Facility* fac, 
                                std::vector<T> value) {
  CheckInput(value);
  value_ = value;
  fac_ = fac;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
FlexibleInput<T>::FlexibleInput(cyclus::Facility* fac,
                                std::vector<T> value, 
                                std::vector<int> time) {
  CheckInput(value, time);
  value_ = value;
  time_ = time;
  fac_ = fac;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
FlexibleInput<T>::~FlexibleInput() {
  delete fac_;
}

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
T FlexibleInput<T>::UpdateValue() {
  // Get current time with t = 0 being the entrance of fac_ in the simulation
  int t = fac_->context()->time() - fac_->enter_time();
  
  if (t >= time_[current_idx_] && t < time_[current_idx_+1]) {
    return value_[current_idx_];
  } else if (t == time_[current_idx_+1]) {
    ++current_idx_;
    return value_[current_idx_];
  } else {
    std::stringstream ss;
    ss << "Agent '" << fac_->prototype()  << "' of spec '" << fac_->spec() 
       << "' with enter_time '" << fac_->enter_time() 
       << "' has passed the invalid timestamp '" << t << "' at time'" 
       << fac_->context()->time() << "' to a FlexibleInput variable.\n";
    
    throw cyclus::ValueError(ss.str());
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
void FlexibleInput<T>::CheckInput_(const std::vector<T>& value) {
  int lifetime = fac_->lifetime();
  if (lifetime == -1) {
    lifetime = fac_->context()->sim_info().duration;
  }

  if (value.size() > lifetime) {
    std::stringstream ss;
    ss << "While initialising agent '" << fac_->prototype()  
       << "' of spec '" << fac_->spec() << "' at time'" 
       << fac_->context()->time() << "' a problem appeared:\n"
       << "The value vector passed to FlexibleInput contains too many "
       << "elements.\n";
    
    throw cyclus::ValueError(ss.str());
 }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <class T>
void FlexibleInput<T>::CheckInput_(const std::vector<T>& value,
                              const std::vector<int>& time) {
  CheckInput_(value);
  if (value.size() != time.size()) {
    std::stringstream ss;
    ss << "While initialising agent '" << fac_->prototype()  
       << "' of spec '" << fac_->spec() << "' at time'" 
       << fac_->context()->time() << "' a problem appeared:\n"
       << "time and value vectors do not have the same size.\n"
       << "size of time vector: " << time.size()
       << ", size of value vector: " << value.size() << "\n";
    
    throw cyclus::ValueError(ss.str());
  }
}

}  // namespace misoenrichment
