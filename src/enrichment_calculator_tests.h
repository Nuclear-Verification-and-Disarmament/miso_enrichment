#ifndef MISOENRICHMENT_SRC_ENRICHMENT_CALCULATOR_TESTS_H_
#define MISOENRICHMENT_SRC_ENRICHMENT_CALCULATOR_TESTS_H_

#include <gtest/gtest.h>

#include <boost/shared_ptr.hpp>

#include "facility_tests.h"
#include "agent_tests.h"
#include "test_context.h"
#include "env.h"
#include "exchange_context.h"

#include "composition.h"
#include "enrichment_calculator.h"

namespace misoenrichment {

class EnrichmentCalculatorTest : public ::testing::Test {
 protected:
  EnrichmentCalculatorTest();

  EnrichmentCalculator e;

  // Values to be calculated by EnrichmentCalculator
  cyclus::Composition::Ptr product_comp, tails_comp;
  cyclus::CompMap product_cm, tails_cm;
  double feed_qty, product_qty, tails_qty, swu_used;
  int n_enriching, n_stripping;

  const cyclus::CompMap expect_product_comp;
  const cyclus::CompMap expect_tails_comp;

  const int expect_n_enriching;
  const int expect_n_stripping;

  const double expect_feed_qty;
  const double expect_product_qty;
  const double expect_tails_qty;
  const double expect_swu_used;

};

}  // namespace misoenrichment

#endif  // MISOENRICHMENT_SRC_ENRICHMENT_CALCULATOR_TESTS_H_
