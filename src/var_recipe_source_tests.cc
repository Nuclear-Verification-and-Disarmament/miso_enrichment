#include "var_recipe_source_tests.h"

#include <gtest/gtest.h>

#include <sstream>

#include "cyc_limits.h"
#include "resource_helpers.h"
#include "test_context.h"

#include "miso_helper.h"

namespace misoenrichment {

void VarRecipeSourceTest::SetUp() {
  simdur = 3;
  src_facility = new VarRecipeSource(tc.get());
  trader = tc.trader();

  InitParameters();
  SetUpVarRecipeSource();
}

void VarRecipeSourceTest::TearDown() {
  delete src_facility;
}

void VarRecipeSourceTest::InitParameters() {
  out_commod = "commod";
  capacity = 42;  // some magic number..

  NucDist nuc_234("uniform", std::vector<double>({5.1e-5, 5.4e-5}));
  NucDist nuc_235("normal", std::vector<double>({7.1e-3, 0.01e-3, 0., 1.}));
  NucDist nuc_238("normalisation", std::vector<double>());
  nuc_map[922340000] = nuc_234;
  nuc_map[922350000] = nuc_235;
  nuc_map[922380000] = nuc_238;
  var_out_recipe = VarOutRecipe("atom", nuc_map);
}

void VarRecipeSourceTest::SetUpVarRecipeSource() {
  src_facility->out_commod = out_commod;
  src_facility->var_out_recipe_ = var_out_recipe;
  src_facility->throughput_times = std::vector<int>({0});
  src_facility->throughput_vals = std::vector<double>({capacity});
  src_facility->current_throughput = capacity;
}

boost::shared_ptr< cyclus::ExchangeContext<cyclus::Material> >
VarRecipeSourceTest::GetContext(int nreqs, std::string commod) {
  using cyclus::Material;
  using cyclus::Request;
  using cyclus::ExchangeContext;
  using test_helpers::get_mat;

  double qty = 3;
  boost::shared_ptr< ExchangeContext<Material> >
      ec(new ExchangeContext<Material>());
  for (int i = 0; i < nreqs; i++) {
    ec->AddRequest(Request<Material>::Create(get_mat(), trader, commod));
  }
  return ec;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(VarRecipeSourceTest, CheckVarOutRecipe_) {
  cyclus::MockSim* fake_sim = new cyclus::MockSim(simdur);
  src_facility = new VarRecipeSource(fake_sim->context());
  SetUpVarRecipeSource();

  EXPECT_NO_THROW(check_var_out_recipe(src_facility));

  var_out_recipe.first = "abc";
  update_var_out_recipe(src_facility);
  EXPECT_THROW(check_var_out_recipe(src_facility), cyclus::ValueError);

  var_out_recipe.first = "mass";
  var_out_recipe.second[922380000].first = "mass";
  update_var_out_recipe(src_facility);
  EXPECT_THROW(check_var_out_recipe(src_facility), cyclus::ValueError);

  var_out_recipe.second[922350000].first = "normalisation";
  var_out_recipe.second[922380000].first = "normalisation";
  update_var_out_recipe(src_facility);
  EXPECT_THROW(check_var_out_recipe(src_facility), cyclus::ValueError);
}

TEST_F(VarRecipeSourceTest, EnterNotify) {
  cyclus::MockSim* fake_sim = new cyclus::MockSim(simdur);
  src_facility = new VarRecipeSource(fake_sim->context());
  SetUpVarRecipeSource();

  EXPECT_NO_THROW(src_facility->EnterNotify());
}

TEST_F(VarRecipeSourceTest, Tick) {
  cyclus::MockSim* fake_sim = new cyclus::MockSim(simdur);
  src_facility = new VarRecipeSource(fake_sim->context());
  SetUpVarRecipeSource();

  // EnterNotify must be called to set up the flexible throughput, else a
  // segfault occurs.
  src_facility->EnterNotify();
  EXPECT_NO_THROW(src_facility->Tick());
}

TEST_F(VarRecipeSourceTest, Tock) {
  cyclus::MockSim* fake_sim = new cyclus::MockSim(simdur);
  src_facility = new VarRecipeSource(fake_sim->context());
  SetUpVarRecipeSource();

  // EnterNotify must be called to set up the flexible throughput, else a
  // segfault occurs.
  src_facility->EnterNotify();
  EXPECT_NO_THROW(src_facility->Tock());
}

TEST_F(VarRecipeSourceTest, AddBids) {
  using cyclus::Bid;
  using cyclus::BidPortfolio;
  using cyclus::CapacityConstraint;
  using cyclus::ExchangeContext;
  using cyclus::Material;

  int nreqs = 5;

  boost::shared_ptr< cyclus::ExchangeContext<Material> >
      ec = GetContext(nreqs, out_commod);

  std::set<BidPortfolio<Material>::Ptr> ports =
      src_facility->GetMatlBids(ec.get()->commod_requests);

  ASSERT_TRUE(ports.size() > 0);
  EXPECT_EQ(ports.size(), 1);

  BidPortfolio<Material>::Ptr port = *ports.begin();
  EXPECT_EQ(port->bidder(), src_facility);
  EXPECT_EQ(port->bids().size(), nreqs);

  const std::set< CapacityConstraint<Material> >& constrs = port->constraints();
  ASSERT_TRUE(constrs.size() > 0);
  EXPECT_EQ(constrs.size(), 1);
  EXPECT_EQ(*constrs.begin(), CapacityConstraint<Material>(capacity));
}

TEST_F(VarRecipeSourceTest, Clone) {
  cyclus::Context* ctx = tc.get();
  misoenrichment::VarRecipeSource* cloned_fac =
    dynamic_cast<misoenrichment::VarRecipeSource*>(src_facility->Clone());

  EXPECT_EQ(outcommod(src_facility),  outcommod(cloned_fac));
  EXPECT_EQ(throughput_times(src_facility), throughput_times(cloned_fac));
  EXPECT_EQ(throughput_vals(src_facility), throughput_vals(cloned_fac));
  EXPECT_EQ(varoutrecipe(src_facility), varoutrecipe(cloned_fac));

  delete cloned_fac;
}

TEST_F(VarRecipeSourceTest, PositionInitialize) {
  std::string config =
      "<out_commod>spent_fuel</out_commod>"
      "<var_out_recipe>"
      "  <mass_or_atom>atom</mass_or_atom>"
      "  <nuclides>"
      "    <item>"
      "      <nuc_id>922340000</nuc_id>"
      "      <rng_properties>"
      "        <distribution>uniform</distribution>"
      "        <parameters>"
      "          <val>0.000051</val>"
      "          <val>0.000054</val>"
      "        </parameters>"
      "      </rng_properties>"
      "    </item>"
      "    <item>"
      "      <nuc_id>922350000</nuc_id>"
      "      <rng_properties>"
      "        <distribution>normal</distribution>"
      "        <parameters>"
      "          <val>0.0071</val>"
      "          <val>0.00001</val>"
      "          <val>0.</val>"
      "          <val>1.</val>"
      "        </parameters>"
      "      </rng_properties>"
      "    </item>"
      "    <item>"
      "      <nuc_id>922380000</nuc_id>"
      "      <rng_properties>"
      "        <distribution>normalisation</distribution>"
      "       <parameters></parameters>"
      "      </rng_properties>"
      "    </item>"
      "  </nuclides>"
      "</var_out_recipe>"
      "<throughput_times><val>0</val></throughput_times>"
      "<throughput_vals><val>10000000</val></throughput_vals>";
  cyclus::MockSim sim(cyclus::AgentSpec(":misoenrichment:VarRecipeSource"), config,
                      simdur);
  int id = sim.Run();

  cyclus::QueryResult qr = sim.db().Query("AgentPosition", NULL);
  EXPECT_EQ(qr.GetVal<double>("Latitude"), 0.0);
  EXPECT_EQ(qr.GetVal<double>("Longitude"), 0.0);
}

TEST_F(VarRecipeSourceTest, Print) {
  EXPECT_NO_THROW(std::string s = src_facility->str());
}

TEST_F(VarRecipeSourceTest, RecordPosition) {
  std::string config =
      "<out_commod>out_commod</out_commod>"
      "<var_out_recipe>"
      "  <mass_or_atom>atom</mass_or_atom>"
      "  <nuclides>"
      "    <item>"
      "      <nuc_id>922340000</nuc_id>"
      "      <rng_properties>"
      "        <distribution>uniform</distribution>"
      "        <parameters>"
      "          <val>0.000051</val>"
      "          <val>0.000054</val>"
      "        </parameters>"
      "      </rng_properties>"
      "    </item>"
      "    <item>"
      "      <nuc_id>922350000</nuc_id>"
      "      <rng_properties>"
      "        <distribution>normal</distribution>"
      "        <parameters>"
      "          <val>0.0071</val>"
      "          <val>0.00001</val>"
      "          <val>0.</val>"
      "          <val>1.</val>"
      "        </parameters>"
      "      </rng_properties>"
      "    </item>"
      "    <item>"
      "      <nuc_id>922380000</nuc_id>"
      "      <rng_properties>"
      "        <distribution>normalisation</distribution>"
      "        <parameters></parameters>"
      "      </rng_properties>"
      "    </item>"
      "  </nuclides>"
      "</var_out_recipe>"
      "<throughput_times><val>0</val></throughput_times>"
      "<throughput_vals><val>10000000</val></throughput_vals>"
      "<latitude>-0.01</latitude>"
      "<longitude>0.01</longitude>";
  int simdur = 3;
  cyclus::MockSim sim(cyclus::AgentSpec(":misoenrichment:VarRecipeSource"), config,
                                        simdur);
  int id = sim.Run();

  cyclus::QueryResult qr = sim.db().Query("AgentPosition", NULL);
  EXPECT_EQ(qr.GetVal<double>("Latitude"), -0.01);
  EXPECT_EQ(qr.GetVal<double>("Longitude"), 0.01);
}

TEST_F(VarRecipeSourceTest, Response) {
  using cyclus::Bid;
  using cyclus::Material;
  using cyclus::Request;
  using cyclus::Trade;
  using test_helpers::get_mat;

  cyclus::MockSim* fake_sim = new cyclus::MockSim(simdur);
  src_facility = new VarRecipeSource(fake_sim->context());
  SetUpVarRecipeSource();
  src_facility->EnterNotify();
  src_facility->Tick();

  std::vector< cyclus::Trade<cyclus::Material> > trades;
  std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                        cyclus::Material::Ptr> > responses;

  // Null response
  EXPECT_NO_THROW(src_facility->GetMatlTrades(trades, responses));
  EXPECT_EQ(responses.size(), 0);

  double qty = capacity / 3;
  Request<Material>* request =
      Request<Material>::Create(get_mat(), trader, out_commod);
  Bid<Material>* bid =
      Bid<Material>::Create(request, get_mat(), src_facility);

  Trade<Material> trade(request, bid, qty);
  trades.push_back(trade);

  // 1 trade
  src_facility->GetMatlTrades(trades, responses);
  EXPECT_EQ(responses.size(), 1);
  EXPECT_EQ(responses[0].second->quantity(), qty);
  cyclus::CompMap first_comp = responses[0].second->comp()->atom();

  // 2 trades, total qty = capacity
  trades.push_back(trade);
  responses.clear();
  EXPECT_NO_THROW(src_facility->GetMatlTrades(trades, responses));
  EXPECT_EQ(responses.size(), 2);
  cyclus::CompMap second_comp = responses[1].second->comp()->atom();
  // No tick occured, hence the compositions should be identical!
  EXPECT_TRUE(misotest::CompareCompMap(first_comp, second_comp));

  // Reset!
  src_facility->Tick();
  responses.clear();
  src_facility->GetMatlTrades(trades, responses);
  cyclus::CompMap third_comp = responses[0].second->comp()->atom();
  // Tick occured, the composition should change.
  EXPECT_FALSE(misotest::CompareCompMap(first_comp, third_comp));

  delete request;
  delete bid;
}

} // namespace misoenrichment

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Do Not Touch! Below section required for connection with Cyclus
cyclus::Agent* VarRecipeSourceConstructor(cyclus::Context* ctx) {
  return new misoenrichment::VarRecipeSource(ctx);
}

// Required to get functionality in cyclus agent unit tests library
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED

INSTANTIATE_TEST_SUITE_P(VarRecipeSource, FacilityTests,
                         ::testing::Values(&VarRecipeSourceConstructor));
INSTANTIATE_TEST_SUITE_P(VarRecipeSource, AgentTests,
                         ::testing::Values(&VarRecipeSourceConstructor));

