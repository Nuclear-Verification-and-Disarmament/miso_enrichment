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

  const double eps_comp = 1e-5;

  EnrichmentCalculator e;

  const cyclus::CompMap expect_product_comp;
  const cyclus::CompMap expect_tails_comp;
  //const int expect_n_enrich;
  //const int expect_n_strip;
  
  const int expect_n_enriching;
  const int expect_n_stripping;

  const double expect_product_qty;
  const double expect_tails_qty;
  const double expect_swu_used;
  
  /*
  void EnrichmentOutput(cyclus::CompMap& product_comp, 
                        cyclus::CompMap& tails_comp, int& n_enrich, 
                        int& n_strip);
  */

};

}  // namespace misoenrichment

#endif  // MISOENRICHMENT_SRC_ENRICHMENT_CALCULATOR_TESTS_H_
