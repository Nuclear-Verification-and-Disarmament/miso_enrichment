#include <gtest/gtest.h>

#include "multi_isotope_enrich.h"

#include "agent_tests.h"
#include "context.h"
#include "facility_tests.h"
#include "pyhooks.h"

using multiisotopeenrichment::MultiIsotopeEnrich;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class MultiIsotopeEnrichTest : public ::testing::Test {
 protected:
  cyclus::TestContext tc;
  MultiIsotopeEnrich* facility;

  virtual void SetUp() {
    cyclus::PyStart();
    facility = new MultiIsotopeEnrich(tc.get());
  }

  virtual void TearDown() {
    delete facility;
    cyclus::PyStop();
  }
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MultiIsotopeEnrichTest, InitialState) {
  // Test things about the initial state of the facility here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MultiIsotopeEnrichTest, Print) {
  EXPECT_NO_THROW(std::string s = facility->str());
  // Test MultiIsotopeEnrich specific aspects of the print method here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MultiIsotopeEnrichTest, Tick) {
  ASSERT_NO_THROW(facility->Tick());
  // Test MultiIsotopeEnrich specific behaviors of the Tick function here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MultiIsotopeEnrichTest, Tock) {
  EXPECT_NO_THROW(facility->Tock());
  // Test MultiIsotopeEnrich specific behaviors of the Tock function here
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Do Not Touch! Below section required for connection with Cyclus
cyclus::Agent* MultiIsotopeEnrichConstructor(cyclus::Context* ctx) {
  return new MultiIsotopeEnrich(ctx);
}
// Required to get functionality in cyclus agent unit tests library
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED
INSTANTIATE_TEST_CASE_P(MultiIsotopeEnrich, FacilityTests,
                        ::testing::Values(&MultiIsotopeEnrichConstructor));
INSTANTIATE_TEST_CASE_P(MultiIsotopeEnrich, AgentTests,
                        ::testing::Values(&MultiIsotopeEnrichConstructor));
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
