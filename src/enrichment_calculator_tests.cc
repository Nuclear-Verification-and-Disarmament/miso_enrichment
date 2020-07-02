#include <cmath>
#include <gtest/gtest.h>

#include "composition.h"
#include "cyc_limits.h"

#include "enrichment_calculator.h"

namespace multiisotopeenrichment {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST(EnrichmentCalculatorTest, ConcentrationDifference) {
  EnrichmentCalculator e;
  
  cyclus::CompMap product_comp;
  product_comp[922350000] = 0.2;
  product_comp[922380000] = 0.8;

  cyclus::CompMap tails_comp;
  tails_comp[922350000] = 0.002;
  tails_comp[922380000] = 0.998;

  e.product_composition = product_comp;
  e.tails_composition = tails_comp;
  e.target_product_assay = 0.1;
  e.target_tails_assay = 0.001;
  
  EXPECT_TRUE(cyclus::AlmostEq(2., 
                               std::pow(e.ConcentrationDifference(), 2.)));
}

}  // namespace multiisotopeenrichment

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED

