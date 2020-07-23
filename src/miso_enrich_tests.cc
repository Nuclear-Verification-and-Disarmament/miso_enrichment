#include <gtest/gtest.h>

#include "miso_enrich_tests.h"

#include "agent_tests.h"
#include "context.h"
#include "env.h"
#include "facility_tests.h"

namespace misoenrichment {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrichTest::SetUp() {
  cyclus::Env::SetNucDataPath();
  cyclus::Context* ctx = tc_.get();
  miso_enrich_facility = new MIsoEnrich(ctx);
  trader = tc_.trader();
  InitParameters();
  SetUpMIsoEnrichment();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrichTest::TearDown() {
  delete miso_enrich_facility;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrichTest::InitParameters() {
  cyclus::Context* ctx = tc_.get();

  cyclus::CompMap cm;
  double feed_assay = 0.711;
  double feed_U234 = 5.5e-3;
  cm[922340000] = feed_U234;
  cm[922350000] = feed_assay;
  cm[922380000] = 100 - feed_assay - feed_U234;
  recipe = cyclus::Composition::CreateFromMass(cm);
  feed_recipe = "feed_recipe";
  ctx->AddRecipe(feed_recipe, recipe);

  feed_commod = "feed_U";
  product_commod = "enriched_U";
  tails_commod = "depleted_U";
  tails_assay = 0.001;
  initial_feed = 10;
  inv_size = 100;
  order_prefs = true;
  max_enrich = 0.9;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrichTest::SetUpMIsoEnrichment() {
  miso_enrich_facility->feed_recipe = feed_recipe;
  miso_enrich_facility->feed_commod = feed_commod;
  miso_enrich_facility->product_commod = product_commod;
  miso_enrich_facility->tails_commod = tails_commod;
  miso_enrich_facility->tails_assay = tails_assay;
  miso_enrich_facility->initial_feed = initial_feed;
  miso_enrich_facility->max_feed_inventory = inv_size;
  miso_enrich_facility->order_prefs = order_prefs;
  miso_enrich_facility->swu_capacity = swu_capacity;
  miso_enrich_facility->max_enrich = max_enrich;
  miso_enrich_facility->latitude = latitude;
  miso_enrich_facility->longitude = longitude;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Material::Ptr MIsoEnrichTest::GetFeedMat(double qty) {
  return cyclus::Material::CreateUntracked(
      qty, tc_.get()->GetRecipe(feed_recipe));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, Request) {
  // Test correctness of quantities and materials requested
  double req_qty = inv_size;
  double add_qty;

  cyclus::Material::Ptr mat = miso_enrich_facility->Request_();
  EXPECT_DOUBLE_EQ(mat->quantity(), req_qty);
  EXPECT_EQ(mat->comp(), tc_.get()->GetRecipe(feed_recipe));

  add_qty = inv_size / 3.;
  ASSERT_NO_THROW(miso_enrich_facility->AddMat_(GetFeedMat(add_qty)));
  req_qty -= add_qty;
  mat = miso_enrich_facility->Request_();
  EXPECT_DOUBLE_EQ(mat->quantity(), req_qty);
  EXPECT_EQ(mat->comp(), tc_.get()->GetRecipe(feed_recipe));

  add_qty = 2 * inv_size / 3.;
  ASSERT_NO_THROW(miso_enrich_facility->AddMat_(GetFeedMat(add_qty)));
  req_qty -= add_qty;
  ASSERT_DOUBLE_EQ(req_qty, 0.);
  mat = miso_enrich_facility->Request_();
  EXPECT_DOUBLE_EQ(mat->quantity(), req_qty);
  EXPECT_EQ(mat->comp(), tc_.get()->GetRecipe(feed_recipe));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, ValidRequest) {
  // Check that U238 is present and that the requested U235 enrichment
  // lies between the tails assays and the maximum feasible product assay
  using cyclus::Composition;
  
  //Composition::Ptr depletedU = misoenrichmenttest::depletedU();
  //cyclus::Material::Ptr
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Predefined tests by the cycstub program. TODO delete them or keep them?
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, Print) {
  EXPECT_NO_THROW(std::string s = miso_enrich_facility->str());
  // Test MIsoEnrich specific aspects of the print method here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, Tick) {
  ASSERT_NO_THROW(miso_enrich_facility->Tick());
  // Test MIsoEnrich specific behaviors of the Tick function here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, Tock) {
  EXPECT_NO_THROW(miso_enrich_facility->Tock());
  // Test MIsoEnrich specific behaviors of the Tock function here
}

}  // namespace misoenrichment

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Do Not Touch! Below section required for connection with Cyclus
cyclus::Agent* MIsoEnrichConstructor(cyclus::Context* ctx) {
  return new misoenrichment::MIsoEnrich(ctx);
}
// Required to get functionality in cyclus agent unit tests library
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED
INSTANTIATE_TEST_CASE_P(MIsoEnrich, FacilityTests,
                        ::testing::Values(&MIsoEnrichConstructor));
INSTANTIATE_TEST_CASE_P(MIsoEnrich, AgentTests,
                        ::testing::Values(&MIsoEnrichConstructor));
