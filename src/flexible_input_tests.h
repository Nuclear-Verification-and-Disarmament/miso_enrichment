#ifndef MISOENRICHMENT_SRC_FLEXIBLE_INPUT_TESTS_H_
#define MISOENRICHMENT_SRC_FLEXIBLE_INPUT_TESTS_H_

#include <gtest/gtest.h>

//#include "flexible_input.h"
#include "flexible_input.cc"

namespace cyclus {
  class Agent;
  class MockSim;
}

namespace misoenrichment {

class FlexibleInputTest : public ::testing::Test {
 protected:
  FlexibleInputTest();
  ~FlexibleInputTest();

  cyclus::MockSim SetUpMockSim();
  cyclus::Agent* parent;
  const int duration = 10;
};

}  // namespace misoenrichment

#endif  // MISOENRICHMENT_SRC_FLEXIBLE_INPUT_TESTS_H_
