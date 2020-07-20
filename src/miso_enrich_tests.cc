#include <gtest/gtest.h>

#include "miso_enrich.h"

#include "agent_tests.h"
#include "context.h"
#include "facility_tests.h"
#include "pyhooks.h"

using misoenrichment::MIsoEnrich;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class MIsoEnrichTest : public ::testing::Test {
 protected:
  cyclus::TestContext tc;
  MIsoEnrich* facility;

  virtual void SetUp() {
    cyclus::PyStart();
    facility = new MIsoEnrich(tc.get());
  }

  virtual void TearDown() {
    delete facility;
    cyclus::PyStop();
  }
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, InitialState) {
  // Test things about the initial state of the facility here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, Print) {
  EXPECT_NO_THROW(std::string s = facility->str());
  // Test MIsoEnrich specific aspects of the print method here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, Tick) {
  ASSERT_NO_THROW(facility->Tick());
  // Test MIsoEnrich specific behaviors of the Tick function here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, Tock) {
  EXPECT_NO_THROW(facility->Tock());
  // Test MIsoEnrich specific behaviors of the Tock function here
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Do Not Touch! Below section required for connection with Cyclus
cyclus::Agent* MIsoEnrichConstructor(cyclus::Context* ctx) {
  return new MIsoEnrich(ctx);
}
// Required to get functionality in cyclus agent unit tests library
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED
INSTANTIATE_TEST_CASE_P(MIsoEnrich, FacilityTests,
                        ::testing::Values(&MIsoEnrichConstructor));
INSTANTIATE_TEST_CASE_P(MIsoEnrich, AgentTests,
                        ::testing::Values(&MIsoEnrichConstructor));
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
