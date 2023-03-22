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
    expect_swu_used(199.17105),
    expect_n_enriching(56),
    expect_n_stripping(14) {
  bool use_downblending = false;
  bool use_integer_stages = true;
  double target_product_assay = 0.9;
  double target_tails_assay = 0.001;
  double gamma = 1.3;
  double target_feed_qty = 100;
  double target_product_qty = 1e299;
  double max_swu = 1e299;

  e = EnrichmentCalculator(compPtr_nat_U(), target_product_assay,
                             target_tails_assay, gamma, target_feed_qty,
                             target_product_qty, max_swu,
                             use_downblending, use_integer_stages);
  e.EnrichmentOutput(product_comp, tails_comp, feed_qty, swu_used,
                     product_qty, tails_qty, n_enriching, n_stripping);
  product_cm = product_comp->atom();
  tails_cm = tails_comp->atom();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, DownblendingOptions) {
  double target_product_assay = 0.9;
  double target_tails_assay = 0.001;
  double gamma = 1.3;
  double target_feed_qty = 100;
  double target_product_qty = 1e299;
  double max_swu = 1e299;

  bool use_downblending = true;
  bool use_integer_stages = true;
  EXPECT_NO_THROW(
    EnrichmentCalculator(compPtr_nat_U(), target_product_assay,
                         target_tails_assay, gamma, target_feed_qty,
                         target_product_qty, max_swu,
                         use_downblending, use_integer_stages)
  );
  use_downblending = false;
  EXPECT_NO_THROW(
    EnrichmentCalculator(compPtr_nat_U(), target_product_assay,
                         target_tails_assay, gamma, target_feed_qty,
                         target_product_qty, max_swu,
                         use_downblending, use_integer_stages)
  );
  use_integer_stages = false;
  EXPECT_NO_THROW(
    EnrichmentCalculator(compPtr_nat_U(), target_product_assay,
                         target_tails_assay, gamma, target_feed_qty,
                         target_product_qty, max_swu,
                         use_downblending, use_integer_stages)
  );
  use_downblending = true;
  EXPECT_THROW(
    EnrichmentCalculator(compPtr_nat_U(), target_product_assay,
                         target_tails_assay, gamma, target_feed_qty,
                         target_product_qty, max_swu,
                         use_downblending, use_integer_stages),
    cyclus::ValueError
  );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, AssignmentOperator) {
  EnrichmentCalculator e2(
    cyclus::Composition::CreateFromAtom(weapons_grade_U()), 0.95, 0.1, 1.1,
    1., 1e299, 1e299, true);
  e2 = e;

  cyclus::Composition::Ptr product_comp2, tails_comp2;
  double feed_qty2, product_qty2, tails_qty2, swu_used2;
  double n_enriching2, n_stripping2;

  e2.EnrichmentOutput(product_comp2, tails_comp2, feed_qty2, swu_used2,
                      product_qty2, tails_qty2, n_enriching2, n_stripping2);

  // This test does not strictly check the correct working of the
  // assignment operator, but it does check the results. It is expected
  // that if the assignment should not have worked correctly then the
  // results would be wrong as well.
  EXPECT_TRUE(misotest::CompareCompMap(product_comp2->atom(),
                                       product_comp->atom()));
  EXPECT_TRUE(misotest::CompareCompMap(tails_comp2->atom(),
                                       tails_comp->atom()));
  EXPECT_DOUBLE_EQ(feed_qty2, feed_qty);
  EXPECT_DOUBLE_EQ(product_qty2, product_qty);
  EXPECT_DOUBLE_EQ(tails_qty2, tails_qty);
  EXPECT_DOUBLE_EQ(swu_used2, swu_used);
  EXPECT_DOUBLE_EQ(n_enriching2, n_enriching);
  EXPECT_DOUBLE_EQ(n_stripping2, n_stripping);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, Concentrations) {
  EXPECT_TRUE(misotest::CompareCompMap(expect_product_comp,
                                        product_cm));
  EXPECT_TRUE(misotest::CompareCompMap(expect_tails_comp, tails_cm));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, Flows) {
  EXPECT_NEAR(expect_feed_qty, feed_qty, kEpsDouble);
  EXPECT_NEAR(expect_product_qty, product_qty, kEpsDouble);
  EXPECT_NEAR(expect_tails_qty, tails_qty, kEpsDouble);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, Swu) {
  EXPECT_NEAR(expect_swu_used, swu_used, kEpsDouble);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, NumberStages) {
  EXPECT_DOUBLE_EQ(expect_n_enriching, n_enriching);
  EXPECT_DOUBLE_EQ(expect_n_stripping, n_stripping);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, Downblending) {
  double target_product_assay = MIsoAssay(weapons_grade_U()) - 0.001;

  cyclus::Composition::Ptr bl_product_comp, bl_tails_comp;
  cyclus::Composition::Ptr bl_product_comp2;
  double bl_feed_qty, bl_product_qty, bl_tails_qty, bl_swu_used;
  double bl_feed_qty2, bl_product_qty2;
  double dummy_double;

  // In this case, the feed is the constraining factor.
  EnrichmentCalculator blender(compPtr_nat_U(), target_product_assay,
                               0.001, 1.3, 100, 1e299, 1e299, true);
  blender.EnrichmentOutput(bl_product_comp, bl_tails_comp, bl_feed_qty,
                           bl_swu_used, bl_product_qty, bl_tails_qty,
                           dummy_double, dummy_double);

  EXPECT_DOUBLE_EQ(target_product_assay, MIsoAtomAssay(bl_product_comp));
  EXPECT_DOUBLE_EQ(expect_feed_qty, bl_feed_qty);
  EXPECT_TRUE(bl_product_qty > expect_product_qty);

  // In this case, the product is the constraining factor.
  // Of course, the results from this and from the previous enrichment
  // should be identical.
  EnrichmentCalculator blender2(compPtr_nat_U(), target_product_assay,
                                0.001, 1.3, 1e299, bl_product_qty,
                                1e299, true);
  blender2.EnrichmentOutput(bl_product_comp2, bl_tails_comp, bl_feed_qty2,
                            bl_swu_used, bl_product_qty2, bl_tails_qty,
                            dummy_double, dummy_double);

  EXPECT_DOUBLE_EQ(target_product_assay, MIsoAtomAssay(bl_product_comp2));
  EXPECT_DOUBLE_EQ(bl_feed_qty, bl_feed_qty2);
  EXPECT_DOUBLE_EQ(bl_product_qty, bl_product_qty2);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(EnrichmentCalculatorTest, NonIntegerStagesNumbers) {
  double target_product_assay = 0.9;
  double target_tails_assay = 0.001;
  double gamma = 1.3;
  double target_feed_qty = 100;
  double target_product_qty = 1e299;
  double max_swu = 1e299;
  bool use_downblending = false;
  bool use_integer_stages = false;

  cyclus::Composition::Ptr pc, tc;
  cyclus::CompMap expected_product_comp, expected_tails_comp;
  expected_product_comp[922340000] = 0.00772031;
  expected_product_comp[922350000] = 0.900002;
  expected_product_comp[922380000] = 0.09227769;

  expected_tails_comp[922340000] = 2.54688424e-06;
  expected_tails_comp[922350000] = 1.00001000e-03;
  expected_tails_comp[922380000] = 9.98997443e-01;

  EnrichmentCalculator e2(compPtr_nat_U(), target_product_assay,
                          target_tails_assay, gamma, target_feed_qty,
                          target_product_qty, max_swu,
                          use_downblending, use_integer_stages);

  e2.EnrichmentOutput(pc, tc, feed_qty, swu_used,
                      product_qty, tails_qty, n_enriching, n_stripping);
  product_cm = pc->atom();
  tails_cm = tc->atom();
  EXPECT_TRUE(misotest::CompareCompMap(product_cm, expected_product_comp));
  EXPECT_TRUE(misotest::CompareCompMap(tails_cm, expected_tails_comp));
}

}  // namespace misoenrichment

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED

