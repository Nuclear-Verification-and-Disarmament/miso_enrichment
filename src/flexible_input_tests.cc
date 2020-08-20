#include "flexible_input_tests.h"

#include <array>
#include <string>
#include <vector>

#include "dynamic_module.h"  // for cyclus::AgentSpec
#include "error.h"
#include "mock_sim.h"
#include "pyhooks.h"

namespace misoenrichment {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
FlexibleInputTest::FlexibleInputTest() {
  cyclus::PyStart();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
FlexibleInputTest::~FlexibleInputTest() {
  cyclus::PyStop();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::MockSim FlexibleInputTest::SetUpMockSim() {
  std::string agent = ":agents:Source";
  std::string config = "<commod>test_commod</commod>"
                       "<recipe_name>test_recipe</recipe_name>"
                       "<capacity>1</capacity>";

  return cyclus::MockSim(cyclus::AgentSpec(agent), config, duration);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(FlexibleInputTest, CheckSingleInput) {
  cyclus::MockSim sim = SetUpMockSim();
  parent = sim.agent;

  std::vector<int> vals(duration);
  for (int i = 0; i < vals.size(); ++i) {
    vals[i] = i * 10;
  }

  // OK
  EXPECT_NO_THROW(FlexibleInput<int> f(parent, vals););
  
  std::vector<int> wrong_vals = vals;
  wrong_vals.push_back(10);
  
  // Not OK, wrong_vals has more entries than simulation duration
  EXPECT_THROW(FlexibleInput<int> ff(parent, wrong_vals);, 
               cyclus::ValueError);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(FlexibleInputTest, CheckDoubleInput) {
  cyclus::MockSim sim = SetUpMockSim();
  parent = sim.agent;

  const int arr_size = 4;
  int time_arr[arr_size] = {0, 4, 5, 8};
  int val_arr[arr_size] = {10, 20, 30, 40};

  std::vector<int> time(time_arr, time_arr+arr_size);
  std::vector<int> vals(val_arr, val_arr+arr_size);

  // OK
  EXPECT_NO_THROW(FlexibleInput<int> f(parent, vals, time););
  
  vals.pop_back();
  // Not OK, vector sizes mismatch
  EXPECT_THROW(FlexibleInput<int> ff(parent, vals, time);,
               cyclus::ValueError);

  time.erase(time.begin());
  // Not OK, time[0] is not 0
  EXPECT_THROW(FlexibleInput<int> fff(parent, vals, time);,
               cyclus::ValueError);
}

}  // namespace misoenrichment

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED
