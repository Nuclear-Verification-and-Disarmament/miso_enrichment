#include "miso_enrich_tests.h"

#include <set>
#include <vector>

#include "agent_tests.h"
#include "context.h"
#include "dynamic_module.h"
#include "env.h"
#include "facility_tests.h"
#include "mock_sim.h"
#include "pyhooks.h"
#include "query_backend.h"

#include "miso_helper.h"

using cyclus::Cond;
using cyclus::QueryResult;

namespace misoenrichment {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrichTest::SetUp() {
  cyclus::PyStart();
  cyclus::Env::SetNucDataPath();

  fake_sim = new cyclus::MockSim(1);
  miso_enrich_facility = new MIsoEnrich(fake_sim->context());

  InitParameters();
  SetUpMIsoEnrichment();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrichTest::TearDown() {
  delete miso_enrich_facility;
  delete fake_sim;

  cyclus::PyStop();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrichTest::InitParameters() {
  cyclus::CompMap cm;
  double feed_assay = 0.711;
  double feed_U234 = 5.5e-3;
  cm[922340000] = feed_U234;
  cm[922350000] = feed_assay;
  cm[922380000] = 100 - feed_assay - feed_U234;
  recipe = cyclus::Composition::CreateFromMass(cm);
  feed_recipe = "feed_recipe";
  fake_sim->AddRecipe(feed_recipe, recipe);

  feed_commod = "feed_U";
  product_commod = "enriched_U";
  tails_commod = "depleted_U";
  tails_assay = 0.002;
  initial_feed = 0.;
  inv_size = 1000;
  order_prefs = true;
  max_enrich = 0.8;
  swu_capacity = 1e299;
  swu_vals = std::vector<double>(1,1);
  swu_times = std::vector<int>(1,0);
  gamma_235 = 1.4;
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
  miso_enrich_facility->current_swu_capacity = swu_capacity;
  miso_enrich_facility->max_enrich = max_enrich;
  miso_enrich_facility->gamma_235 = gamma_235;
  miso_enrich_facility->latitude = latitude;
  miso_enrich_facility->longitude = longitude;
  miso_enrich_facility->swu_capacity_vals = swu_vals;
  miso_enrich_facility->swu_capacity_times = swu_times;

  miso_enrich_facility->EnterNotify();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Material::Ptr MIsoEnrichTest::GetFeedMat(double qty) {
  return cyclus::Material::CreateUntracked(qty, recipe);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, BidPrefs) {
  // Test that the facility does adjust preferences for different bids
  using cyclus::Composition;
  using cyclus::Material;

  cyclus::CompMap cm;
  cm[922340000] =  0.0055;
  cm[922350000] =  0.7110;
  cm[922380000] = 99.2835;
  cyclus::compmath::Normalize(&cm);
  Composition::Ptr feed_1 = Composition::CreateFromMass(cm);
  cm[922350000] =  1.7110;
  cm[922380000] = 98.2835;
  cyclus::compmath::Normalize(&cm);
  Composition::Ptr feed_2 = Composition::CreateFromMass(cm);

  std::string config =
    "   <feed_commod>feed_U</feed_commod> "
    "   <feed_recipe>feed1</feed_recipe> "
    "   <product_commod>enriched_U</product_commod> "
    "   <tails_commod>depleted_U</tails_commod> "
    "   <tails_assay>0.002</tails_assay> "
    "   <max_feed_inventory>1</max_feed_inventory> "
    "   <order_prefs>1</order_prefs>"
    "   <swu_capacity_times><val>0</val></swu_capacity_times> "
    "   <swu_capacity_vals><val>10000</val></swu_capacity_vals> ";

  int simdur = 1;
  cyclus::MockSim sim(cyclus::AgentSpec(":misoenrichment:MIsoEnrich"),
                      config, simdur);
  sim.AddRecipe("feed1", feed_1);
  sim.AddRecipe("feed2", feed_2);

  sim.AddSource("feed_U").recipe("feed1")
                         .capacity(1)
                         .Finalize();
  sim.AddSource("feed_U").recipe("feed2")
                         .capacity(1)
                         .Finalize();
  int id = sim.Run();

  std::vector<Cond> conds;
  conds.push_back(Cond("Commodity", "==", std::string("feed_U")));
  QueryResult qr = sim.db().Query("Transactions", &conds);

  EXPECT_EQ(1, qr.rows.size());

  Material::Ptr mat = sim.GetMaterial(qr.GetVal<int>("ResourceId"));
  cyclus::CompMap actual =  mat->comp()->mass();
  cyclus::compmath::Normalize(&actual);
  cyclus::compmath::Normalize(&cm);
  misotest::CompareCompMap(actual, cm);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, FeedConstraint) {
  // Check that the feed constraint is evaluated correctly. Only 100 kg of
  // feed are at the disposal but the sink requests infinity kg of product.
  std::string config =
    "   <feed_commod>feed_U</feed_commod> "
    "   <feed_recipe>feed_recipe</feed_recipe> "
    "   <initial_feed>100</initial_feed> "
    "   <product_commod>enriched_U</product_commod> "
    "   <tails_commod>depleted_U</tails_commod> "
    "   <tails_assay>0.002</tails_assay> "
    "   <swu_capacity_times><val>0</val></swu_capacity_times> "
    "   <swu_capacity_vals><val>10000</val></swu_capacity_vals> "
    "   <use_downblending>0</use_downblending> ";

  int simdur = 1;
  cyclus::MockSim sim(cyclus::AgentSpec(":misoenrichment:MIsoEnrich"),
                      config, simdur);
  sim.AddRecipe(feed_recipe, recipe);
  sim.AddRecipe("enriched_U_recipe", misotest::comp_weapongradeU());
  sim.AddSink("enriched_U").recipe("enriched_U_recipe")
                           .Finalize();
  int id = sim.Run();

  std::vector<Cond> conds;
  conds.push_back(Cond("Commodity", "==", std::string("enriched_U")));
  QueryResult qr = sim.db().Query("Transactions", &conds);
  Material::Ptr m = sim.GetMaterial(qr.GetVal<int>("ResourceId"));

  EXPECT_EQ(qr.rows.size(), 1);
  EXPECT_NEAR(m->quantity(),0.5754, 1e-4);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, GetMatlBids) {
  // Test the bidding. At first no bids are expected because there are
  // neither feed nor tails present. Then, feed is added and one bid for
  // the product is expected. Finally, an enrichment is performed and two
  // bids are expected (one for the product, one for the tails).
  using cyclus::Material;

  cyclus::Request<Material> *req_prod, *req_tails;
  cyclus::CommodMap<Material>::type out_requests;
  std::set<cyclus::BidPortfolio<Material>::Ptr> bids;

  cyclus::CompMap cm;
  cm[922340000] = 5.5e-3;
  cm[922350000] = 50;
  cm[922380000] = 49.9945;
  Material::Ptr product = Material::CreateUntracked(
      1, cyclus::Composition::CreateFromMass(cm));
  // Clearing the CompMap is technically not needed but may be clearer.
  cm.clear();
  cm[922350000] = 0.3;
  cm[922380000] = 99.7;
  Material::Ptr tails = Material::CreateUntracked(
      1, cyclus::Composition::CreateFromMass(cm));
  req_prod = cyclus::Request<Material>::Create(product, miso_enrich_facility,
                                               product_commod);
  req_tails = cyclus::Request<Material>::Create(tails, miso_enrich_facility,
                                                tails_commod);

  out_requests[req_prod->commodity()].push_back(req_prod);
  out_requests[req_tails->commodity()].push_back(req_tails);

  ASSERT_NO_THROW(bids = miso_enrich_facility->GetMatlBids(out_requests));
  EXPECT_EQ(0, bids.size());

  DoAddMat(GetFeedMat(1000));
  bids = miso_enrich_facility->GetMatlBids(out_requests);
  EXPECT_EQ(1, bids.size());

  DoEnrich(product, product->quantity());
  bids = miso_enrich_facility->GetMatlBids(out_requests);
  EXPECT_EQ(2, bids.size());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, GetMatlRequests) {
  // Test the size of the RequestPortfolio and that the request asks for
  // the right quantity and commodity.
  std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr> req_ports;
  cyclus::RequestPortfolio<cyclus::Material>::Ptr req_port;
  cyclus::Request<cyclus::Material>* request;
  ASSERT_NO_THROW(req_ports = miso_enrich_facility->GetMatlRequests());

  req_port = *req_ports.begin();
  request = req_port->requests()[0];
  EXPECT_EQ(1, req_ports.size());
  EXPECT_EQ(miso_enrich_facility, req_port->requester());
  EXPECT_EQ(feed_commod, request->commodity());
  EXPECT_DOUBLE_EQ(inv_size - initial_feed, req_port->qty());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, NoBidPrefs) {
  // Test that the facility does adjust preferences for different bids
  using cyclus::Composition;
  using cyclus::Material;

  cyclus::CompMap cm;
  cm[922340000] =  0.0055;
  cm[922350000] =  0.7110;
  cm[922380000] = 99.2835;
  Composition::Ptr feed_1 = Composition::CreateFromMass(cm);
  cm[922350000] =  1.7110;
  cm[922380000] = 98.2835;
  Composition::Ptr feed_2 = Composition::CreateFromMass(cm);

  std::string config =
    "   <feed_commod>feed_U</feed_commod> "
    "   <feed_recipe>feed1</feed_recipe> "
    "   <product_commod>enriched_U</product_commod> "
    "   <tails_commod>depleted_U</tails_commod> "
    "   <tails_assay>0.002</tails_assay> "
    "   <max_feed_inventory>2</max_feed_inventory> "
    "   <order_prefs>1</order_prefs>"
    "   <swu_capacity_times><val>0</val></swu_capacity_times> "
    "   <swu_capacity_vals><val>10000</val></swu_capacity_vals> "
    "   <use_downblending>0</use_downblending> ";
  int simdur = 1;
  cyclus::MockSim sim(cyclus::AgentSpec(":misoenrichment:MIsoEnrich"),
                      config, simdur);
  sim.AddRecipe("feed1", feed_1);
  sim.AddRecipe("feed2", feed_2);

  sim.AddSource("feed_U").recipe("feed1")
                         .capacity(1)
                         .Finalize();
  sim.AddSource("feed_U").recipe("feed2")
                         .capacity(1)
                         .Finalize();
  int id = sim.Run();

  std::vector<Cond> conds;
  conds.push_back(Cond("Commodity", "==", std::string("feed_U")));
  QueryResult qr = sim.db().Query("Transactions", &conds);

  EXPECT_EQ(2, qr.rows.size());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, Print) {
  EXPECT_NO_THROW(std::string s = miso_enrich_facility->str());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, Request) {
  // Test correctness of quantities and materials requested
  double req_qty = inv_size - initial_feed;
  double add_qty;

  // Only initially added material in inventory
  cyclus::Material::Ptr mat = DoRequest();
  EXPECT_DOUBLE_EQ(mat->quantity(), req_qty);
  EXPECT_EQ(mat->comp(), recipe);

  // Full inventory
  add_qty = req_qty;
  ASSERT_NO_THROW(DoAddMat(GetFeedMat(add_qty)));
  req_qty -= add_qty;
  ASSERT_DOUBLE_EQ(req_qty, 0.);
  mat = DoRequest();
  EXPECT_DOUBLE_EQ(mat->quantity(), req_qty);
  EXPECT_EQ(mat->comp(), recipe);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, RequestSim) {
  // Check that exactly as much material as needed is requested in exactly
  // one request.

  std::string config =
    "   <feed_commod>feed_U</feed_commod> "
    "   <feed_recipe>feed_recipe</feed_recipe> "
    "   <product_commod>enriched_U</product_commod> "
    "   <tails_commod>depleted_U</tails_commod> "
    "   <tails_assay>0.002</tails_assay> "
    "   <max_feed_inventory>100</max_feed_inventory> "
    "   <max_enrich>0.8</max_enrich> "
    "   <swu_capacity_times><val>0</val></swu_capacity_times> "
    "   <swu_capacity_vals><val>10000</val></swu_capacity_vals> "
    "   <use_downblending>0</use_downblending> ";
  int simdur = 1;
  cyclus::MockSim sim(cyclus::AgentSpec(":misoenrichment:MIsoEnrich"),
                      config, simdur);

  sim.AddRecipe(feed_recipe, recipe);
  sim.AddSource("feed_U").recipe("feed_recipe").Finalize();

  int id = sim.Run();

  std::vector<Cond> conds;
  conds.push_back(Cond("Commodity", "==", std::string("feed_U")));
  QueryResult qr = sim.db().Query("Transactions", &conds);
  Material::Ptr m = sim.GetMaterial(qr.GetVal<int>("ResourceId"));

  EXPECT_EQ(qr.rows.size(), 1);
  EXPECT_NEAR(m->quantity(), 100, 1e-10);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, TailsTrade) {
  std::string config =
    "   <feed_commod>feed_U</feed_commod> "
    "   <feed_recipe>feed_recipe</feed_recipe> "
    "   <product_commod>enriched_U</product_commod> "
    "   <tails_commod>depleted_U</tails_commod> "
    "   <tails_assay>0.002</tails_assay> "
    "   <initial_feed>100</initial_feed> "
    "   <swu_capacity_times><val>0</val></swu_capacity_times> "
    "   <swu_capacity_vals><val>10000</val></swu_capacity_vals> "
    "   <use_downblending>0</use_downblending> ";
  int simdur = 2;
  cyclus::MockSim sim(cyclus::AgentSpec(":misoenrichment:MIsoEnrich"),
                      config, simdur);
  sim.AddRecipe(feed_recipe, recipe);
  sim.AddRecipe("depleted_U_recipe", misotest::comp_depletedU());
  sim.AddRecipe("enriched_U_recipe", misotest::comp_weapongradeU());

  sim.AddSink("enriched_U").recipe("enriched_U_recipe")
                           .Finalize();
  sim.AddSink("depleted_U").recipe("depleted_U_recipe")
                           .Finalize();

  int id = sim.Run();

  std::vector<Cond> conds;
  conds.push_back(Cond("Commodity", "==", std::string("depleted_U")));
  QueryResult qr = sim.db().Query("Transactions", &conds);
  cyclus::Material::Ptr mat = sim.GetMaterial(
      qr.GetVal<int>("ResourceId"));

  EXPECT_EQ(1, qr.rows.size());
  EXPECT_NEAR(99.4246, mat->quantity(), 1e-4);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, Tick) {
  EXPECT_NO_THROW(miso_enrich_facility->Tick());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, Tock) {
  EXPECT_NO_THROW(miso_enrich_facility->Tock());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(MIsoEnrichTest, ValidRequest) {
  // Check that U238 is present and that the requested U235 enrichment
  // lies between the tails assays and the maximum feasible product assay
  using cyclus::Composition;
  using cyclus::Material;
  using namespace misotest;  // for the predefined compositions

  cyclus::CompMap cm;
  cm[922340000] = 0.5;
  cm[922350000] = 0.5;
  Composition::Ptr comp_no_U238 = Composition::CreateFromMass(cm);

  Material::Ptr naturalU = Material::CreateUntracked(1, recipe);
  Material::Ptr depletedU = Material::CreateUntracked(1, comp_depletedU());
  Material::Ptr no_U238 = Material::CreateUntracked(1, comp_no_U238);
  Material::Ptr weapongradeU = Material::CreateUntracked(
      1, comp_weapongradeU());

  EXPECT_TRUE(DoValidReq(naturalU));
  EXPECT_FALSE(DoValidReq(depletedU));
  EXPECT_FALSE(DoValidReq(no_U238));
  EXPECT_FALSE(DoValidReq(weapongradeU));
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
