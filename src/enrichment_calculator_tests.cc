#include "enrichment_calculator_tests.h"

#include <cmath>
#include <iostream>

#include <gtest/gtest.h>

#include "cyc_limits.h"
#include "comp_math.h"

#include "miso_helper.h"

namespace misoenrichment {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Composition::Ptr compPtr_nat_U() {
  cyclus::CompMap comp;
  comp[922340000] = 5.5e-3;
  comp[922350000] = 0.711;
  comp[922380000] = 100 - comp[922340000] - comp[922350000];
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

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
EnrichmentCalculatorTest::EnrichmentCalculatorTest() :
    expect_product_comp(weapons_grade_U()), 
    expect_tails_comp(depleted_U()), 
    expect_feed_qty(100),
    expect_product_qty(0.67202), 
    expect_tails_qty(99.32798),
    expect_swu_used(199.171),
    expect_n_enriching(56),
    expect_n_stripping(14) {
  double target_product_assay = 0.9;
  double target_tails_assay = 0.001;
  double gamma = 1.3;
  double target_feed_qty = 100;
  double target_product_qty = 1e299;
  double max_swu = 1e299;

  e = EnrichmentCalculator(compPtr_nat_U(), target_product_assay, 
                             target_tails_assay, gamma, target_feed_qty,
                             target_product_qty, max_swu);
  e.EnrichmentOutput(product_comp, tails_comp, feed_qty, swu_used, 
                     product_qty, tails_qty, n_enriching, n_stripping);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, AssignmentOperator) {
  EnrichmentCalculator e2(
    cyclus::Composition::CreateFromAtom(weapons_grade_U()), 0.95, 0.1, 1.1,
    1.);
  e2 = e;

  cyclus::CompMap product_comp2, tails_comp2;
  double feed_qty2, product_qty2, tails_qty2, swu_used2;
  int n_enriching2, n_stripping2;
  
  e2.EnrichmentOutput(product_comp2, tails_comp2, feed_qty2, swu_used2, 
                      product_qty2, tails_qty2, n_enriching2, n_stripping2);
  
  // This test does not strictly check the correct working of the 
  // assignment operator, but it does check the results. It is expected
  // that if the assignment should not have worked correctly then the 
  // results would be wrong as well.
  EXPECT_TRUE(misotest::CompareCompMap(product_comp2, product_comp));
  EXPECT_TRUE(misotest::CompareCompMap(tails_comp2, tails_comp));
  EXPECT_DOUBLE_EQ(feed_qty2, feed_qty);
  EXPECT_DOUBLE_EQ(product_qty2, product_qty);
  EXPECT_DOUBLE_EQ(tails_qty2, tails_qty);
  EXPECT_DOUBLE_EQ(swu_used2, swu_used);
  EXPECT_EQ(n_enriching2, n_enriching);
  EXPECT_EQ(n_stripping2, n_stripping);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, Concentrations) {
  EXPECT_TRUE(misotest::CompareCompMap(expect_product_comp, 
                                        product_comp));
  EXPECT_TRUE(misotest::CompareCompMap(expect_tails_comp, tails_comp));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, Flows) {
  const double abs_tol = 1e-5;  // maximum absolute difference
  EXPECT_NEAR(expect_feed_qty, feed_qty, abs_tol);
  EXPECT_NEAR(expect_product_qty, product_qty, abs_tol);
  EXPECT_NEAR(expect_tails_qty, tails_qty, abs_tol);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, Swu) {
  const double abs_tol = 1e-5;  // maximum absolute difference
  EXPECT_NEAR(expect_swu_used, swu_used, abs_tol);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, NumberStages) {
  EXPECT_EQ(expect_n_enriching, n_enriching);
  EXPECT_EQ(expect_n_stripping, n_stripping);
}

}  // namespace misoenrichment

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED

