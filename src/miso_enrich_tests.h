#ifndef MISOENRICHMENT_SRC_MISO_ENRICH_TESTS_H_
#define MISOENRICHMENT_SRC_MISO_ENRICH_TESTS_H_

#include <string>

#include <gtest/gtest.h>

#include "composition.h"
#include "material.h"

#include "miso_enrich.h"

namespace misoenrichment {

class MIsoEnrichTest : public ::testing::Test {
 protected:
  MIsoEnrich* miso_enrich_facility;
  
  cyclus::MockSim* fake_sim;
  cyclus::Composition::Ptr recipe;
  
  bool order_prefs;
  double gamma_235, initial_feed, inv_size, latitude, longitude,
         max_enrich, swu_capacity, tails_assay;

  std::string feed_commod, product_commod, tails_commod, feed_recipe;
  std::vector<double> swu_vals;
  std::vector<int> swu_times;
  
  // Functions to initialise and end each test
  void SetUp();
  void TearDown();
  void InitParameters();
  void SetUpMIsoEnrichment();
  
  // Helper functions
  cyclus::Material::Ptr GetFeedMat(double qty);

  // The Do* functions are a hack: MIsoEnrichTest is a friend class to 
  // MIsoEnrich and it can therefore access private functions. However,
  // all of the tests are sub-classes of the fixture and they cannot access
  // the private functions, hence we need to use the Do* functions as an
  // intermediary. See, e.g., 4th bullet point here:
  // https://github.com/google/googletest/blob/master/googletest/docs/advanced.md#testing-private-code
  inline bool DoValidReq(const cyclus::Material::Ptr mat) {
    return miso_enrich_facility->ValidReq_(mat);
  }
  inline cyclus::Material::Ptr DoRequest() {
    return miso_enrich_facility->Request_();
  }
  inline void DoAddMat(cyclus::Material::Ptr mat) {
    miso_enrich_facility->AddMat_(mat);
  }
  inline void DoAddFeedMat(cyclus::Material::Ptr mat) { 
    miso_enrich_facility->AddFeedMat_(mat);
  }
  inline void DoEnrich(cyclus::Material::Ptr mat, double qty) {
    miso_enrich_facility->Enrich_(mat, qty);
  }
};

}  // namespace misoenrichment

#endif  // MISOENRICHMENT_SRC_MISO_ENRICH_TESTS_H_
