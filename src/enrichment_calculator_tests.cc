#include "enrichment_calculator_tests.h"

#include <cmath>
#include <iostream>
#include <gtest/gtest.h>

#include "cyc_limits.h"
#include "comp_math.h"

#include "multi_isotope_helper.h"

namespace multiisotopeenrichment {

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
    expect_product_qty(0.67202), 
    expect_tails_qty(99.32798),
    expect_swu_used(199.171),
    expect_n_enriching(56),
    expect_n_stripping(14) {
  double target_product_assay = 0.9;
  double target_tails_assay = 0.001;
  double gamma = 1.3;
  double feed_qty = 100;

  e = EnrichmentCalculator(compPtr_nat_U(), target_product_assay, 
                             target_tails_assay, gamma, feed_qty);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bool EnrichmentCalculatorTest::expect_true_compmap(cyclus::CompMap cm1,
                                                   cyclus::CompMap cm2) {
  std::vector<int> isotopes;
  IsotopesNucID(isotopes);
  std::vector<int>::iterator it;
  for (it = isotopes.begin(); it != isotopes.end(); it++) {
    cm1[*it] += 1e-299;
    cm2[*it] += 1e-299;
  }  // The code snippet above has been added to ensure that the all of the
     // uranium keys are present in both compmaps, else the comparison 
     // fails.

  bool result = cyclus::compmath::AlmostEq(cm1, cm2, eps_comp); 
  EXPECT_TRUE(result);
  if (!result) {
    std::cout << "Value of: cm1\n"
              << "Actual:\n";
    cyclus::CompMap::iterator it;
    for (it = cm1.begin(); it!= cm1.end(); it++) {
      std::cout << "          " << it->first << ": " << it->second << "\n";
    }
    std::cout << "Expected: \n";
    for (it = cm2.begin(); it!= cm2.end(); it++) {
      std::cout << "          " << it->first << ": " << it->second << "\n";
    }
    std::cout << "\n";
  }
  return result;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, AssignmentOperator) {
  EnrichmentCalculator e2(
    cyclus::Composition::CreateFromAtom(weapons_grade_U()), 0.95, 0.1, 1.1,
    1.);
  e2 = e;
  e2.BuildMatchedAbundanceRatioCascade();

  cyclus::CompMap product_comp1, product_comp2;
  cyclus::CompMap tails_comp1, tails_comp2;
  double feed_qty1, feed_qty2;
  double product_qty1, product_qty2;
  double tails_qty1, tails_qty2;
  double max_swu1, max_swu2;
  int n_enrich1, n_enrich2;
  int n_strip1, n_strip2;
  
  e.EnrichmentOutput(product_comp1, tails_comp1, feed_qty1, product_qty1,
                     tails_qty1, max_swu1, n_enrich1, n_strip1);
  e2.EnrichmentOutput(product_comp2, tails_comp2, feed_qty2, product_qty2,
                      tails_qty2, max_swu2, n_enrich2, n_strip2);
  
  // This test does not strictly check the correct working of the 
  // assignment operator, but it does check the results. It is expected
  // that if the assignment should not have worked correctly then the 
  // results would be wrong as well.
  expect_true_compmap(product_comp2, product_comp1);
  expect_true_compmap(tails_comp2, tails_comp1);
  EXPECT_DOUBLE_EQ(feed_qty2, feed_qty1);
  EXPECT_DOUBLE_EQ(product_qty2, product_qty1);
  EXPECT_DOUBLE_EQ(tails_qty2, tails_qty1);
  EXPECT_DOUBLE_EQ(max_swu2, max_swu1);
  EXPECT_EQ(n_enrich2, n_enrich1);
  EXPECT_EQ(n_strip2, n_strip1);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, Concentrations) {
  double dummy_double;
  int dummy_int;
  cyclus::CompMap product_comp, tails_comp;
  e.EnrichmentOutput(product_comp, tails_comp, dummy_double, dummy_double,
                     dummy_double, dummy_double, dummy_int, dummy_int);

  expect_true_compmap(expect_product_comp, product_comp);
  expect_true_compmap(expect_tails_comp, tails_comp);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, Flows) {
  cyclus::CompMap dummy_compmap;
  double dummy_double, feed_used, product_qty, tails_qty, swu_used;
  int dummy_int;

  e.EnrichmentOutput(dummy_compmap, dummy_compmap, feed_used, swu_used,
                     product_qty, tails_qty, dummy_int, dummy_int);
  // TODO implement feed_used qty
  EXPECT_DOUBLE_EQ(expect_product_qty, product_qty);
  EXPECT_DOUBLE_EQ(expect_tails_qty, tails_qty);
  EXPECT_DOUBLE_EQ(expect_swu_used, swu_used);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, NumberStages) {
  cyclus::CompMap dummy_compmap;
  double dummy_double;
  int n_enriching, n_stripping;

  e.EnrichmentOutput(dummy_compmap, dummy_compmap, dummy_double,
                     dummy_double, dummy_double, dummy_double, 
                     n_enriching, n_stripping);
  EXPECT_EQ(expect_n_enriching, n_enriching);
  EXPECT_EQ(expect_n_stripping, n_stripping);
}


  /*
  void EnrichmentOutput(
      cyclus::CompMap& product_comp, cyclus::CompMap& tails_comp, 
      double& feed_used, double& swu_used, double& product_produced, 
 */


/* --------------------------------------------------------------------

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
   -------------------------------------------------------------------- */

}  // namespace multiisotopeenrichment

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED

