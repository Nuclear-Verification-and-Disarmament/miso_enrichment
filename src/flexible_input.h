#ifndef MISOENRICHMENT_SRC_FLEXIBLE_INPUT_H_
#define MISOENRICHMENT_SRC_FLEXIBLE_INPUT_H_

#include <vector>

#include "facility.h"

/*
// Forward declaration to be able to access lifetime of the facility and
// the current simulation time
namespace cyclus {

class Facility;

}  // namespace cyclus
*/

namespace misoenrichment {

template <class T>
class FlexibleInput {
 public:
  FlexibleInput(cyclus::Facility* fac, std::vector<T> value);
  FlexibleInput(cyclus::Facility* fac, std::vector<T> value, 
                std::vector<int> time);
  ~FlexibleInput();

  T UpdateValue();

 private:
  void CheckInput_(const std::vector<T>& value);
  void CheckInput_(const std::vector<T>& value, 
                   const std::vector<int>& time);

  std::vector<T> value_;
  std::vector<int> time_;
  
  int current_idx_ = 0;

  cyclus::Facility* fac_;
};

}  // namespace misoenrichment

#endif  // MISOENRICHMENT_SRC_FLEXIBLE_INPUT_H_
