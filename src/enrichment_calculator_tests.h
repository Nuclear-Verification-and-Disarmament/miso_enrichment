#ifndef MULTIISOTOPEENRICHMENT_SRC_ENRICHMENT_CALCULATOR_TESTS_H_
#define MULTIISOTOPEENRICHMENT_SRC_ENRICHMENT_CALCULATOR_TESTS_H_

#include <gtest/gtest.h>

#include "composition.h"
#include "enrichment_calculator.h"

namespace multiisotopeenrichment {

class EnrichmentCalculatorTest : public ::testing::Test {
 protected:
  EnrichmentCalculatorTest();

  EnrichmentCalculator e;

  const cyclus::Composition::Ptr feed_comp;
  const cyclus::CompMap expect_product_comp;
  const cyclus::CompMap expect_tails_comp;
  
  double target_product_assay;
  double target_tails_assay;
  
  double feed_qty;
  double expect_product_qty;
  double expect_tails_qty;
  double gamma;

};

}  // namespace multiisotopeenrichment

#endif  // MULTIISOTOPEENRICHMENT_SRC_ENRICHMENT_CALCULATOR_TESTS_H_
