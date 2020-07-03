#include "enrichment_calculator_tests.h"

#include <cmath>
#include <gtest/gtest.h>

#include "cyc_limits.h"

namespace multiisotopeenrichment {

cyclus::Composition::Ptr compPtr_nat_U() {
  cyclus::CompMap comp;
  comp[922340000] = 5.5e-3;
  comp[922350000] = 0.711;
  comp[922380000] = 100 - comp[92234000] - comp[922350000];
  cyclus::Composition::Ptr compPtr = cyclus::Composition::CreateFromAtom(
      comp);
  return compPtr;
}

cyclus::CompMap weapons_grade_U() {
  cyclus::CompMap comp;
  comp[922340000] = 0.00780791;
  comp[922350000] = 0.91020719;
  comp[922380000] = 0.08198490;
  return comp;
}

cyclus::CompMap depleted_U() {
  cyclus::CompMap comp;
  comp[922340000] = 0.0000025465;
  comp[922350000] = 0.0009999580;
  comp[922380000] = 0.9989974955;
  return comp;
}

EnrichmentCalculatorTest::EnrichmentCalculatorTest() :
    feed_comp(compPtr_nat_U()), 
    expect_product_comp(weapons_grade_U()), 
    expect_tails_comp(depleted_U()), 
    target_product_assay(0.9), target_tails_assay(0.001),
    feed_qty(100), 
    expect_product_qty(0.67202), 
    expect_tails_qty(99.32798),
    gamma(1.3),
    e(EnrichmentCalculator(feed_comp, target_product_assay, 
                           target_tails_assay, gamma, 
                           feed_qty=this->feed_qty))
  {;}


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

