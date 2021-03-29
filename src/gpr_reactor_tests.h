#ifndef MISOENRICHMENT_SRC_GPR_REACTOR_TESTS_H_
#define MISOENRICHMENT_SRC_GPR_REACTOR_TESTS_H_

#include <gtest/gtest.h>

#include "agent_tests.h"
#include "gpr_reactor.h"

namespace misoenrichment {

class GprReactorTest : public ::testing::Test {
 protected:
  cyclus::TestContext tc;
  GprReactor* facility;

  virtual void SetUp();
  virtual void TearDown();
};

}  // namespace misoenrichment

#endif  // MISOENRICHMENT_SRC_GPR_REACTOR_TESTS_H_
