#include "gpr_reactor_tests.h"

#include "agent_tests.h"
#include "context.h"
#include "facility_tests.h"
#include "pyhooks.h"

namespace misoenrichment {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void GprReactorTest::SetUp() {
    cyclus::PyStart();
    facility = new GprReactor(tc.get());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void GprReactorTest::TearDown() {
    delete facility;
    cyclus::PyStop();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(GprReactorTest, InitialState) {
  // Test things about the initial state of the facility here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(GprReactorTest, Print) {
  EXPECT_NO_THROW(std::string s = facility->str());
  // Test GprReactor specific aspects of the print method here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(GprReactorTest, Tick) {
  ASSERT_NO_THROW(facility->Tick());
  // Test GprReactor specific behaviors of the Tick function here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(GprReactorTest, Tock) {
  EXPECT_NO_THROW(facility->Tock());
  // Test GprReactor specific behaviors of the Tock function here
}

}  // namespace misoenrichment

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Do Not Touch! Below section required for connection with Cyclus
cyclus::Agent* GprReactorConstructor(cyclus::Context* ctx) {
  return new misoenrichment::GprReactor(ctx);
}
// Required to get functionality in cyclus agent unit tests library
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED
INSTANTIATE_TEST_CASE_P(GprReactor, FacilityTests,
                        ::testing::Values(&GprReactorConstructor));
INSTANTIATE_TEST_CASE_P(GprReactor, AgentTests,
                        ::testing::Values(&GprReactorConstructor));
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
