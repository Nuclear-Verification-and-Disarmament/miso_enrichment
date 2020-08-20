#ifndef MISOENRICHMENT_SRC_FLEXIBLE_INPUT_H_
#define MISOENRICHMENT_SRC_FLEXIBLE_INPUT_H_

#include <vector>

#include "agent.h"

namespace misoenrichment {

template <class T>
class FlexibleInput {
 public:
  FlexibleInput();
  FlexibleInput(cyclus::Agent* parent, std::vector<T> value);
  FlexibleInput(cyclus::Agent* parent, std::vector<T> value, 
                std::vector<int> time);

  T UpdateValue();

 private:
  void CheckInput_(const std::vector<T>& value);
  void CheckInput_(const std::vector<T>& value, 
                   const std::vector<int>& time);

  std::vector<T> value_;
  std::vector<int> time_;
  
  int current_idx_ = 0;

  cyclus::Agent* parent_;
};

}  // namespace misoenrichment

#endif  // MISOENRICHMENT_SRC_FLEXIBLE_INPUT_H_
