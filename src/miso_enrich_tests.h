#ifndef MISOENRICHMENT_SRC_MISO_ENRICH_TESTS_H_
#define MISOENRICHMENT_SRC_MISO_ENRICH_TESTS_H_

#include <string>

#include <gtest/gtest.h>

#include "composition.h"
#include "material.h"
#include "test_context.h"
#include "test_agents/test_facility.h"

#include "miso_enrich.h"

namespace misoenrichment {

class MIsoEnrichTest : public ::testing::Test {
 protected:
  MIsoEnrich* miso_enrich_facility;

  cyclus::TestContext tc_;
  cyclus::Composition::Ptr recipe;
  
  double initial_feed, inv_size, latitude, longitude, max_enrich,
         swu_capacity, tails_assay;

  bool order_prefs;

  std::string feed_commod, product_commod, tails_commod, feed_recipe;
  
  TestFacility* trader;

  void SetUp();
  void TearDown();
  void InitParameters();
  void SetUpMIsoEnrichment();
  
  cyclus::Material::Ptr GetFeedMat(double qty);

};

}  // namespace misoenrichment

#endif  // MISOENRICHMENT_SRC_MISO_ENRICH_TESTS_H_
