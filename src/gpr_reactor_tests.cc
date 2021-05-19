#include "gpr_reactor_tests.h"

#include <cstdio>

#include "agent_tests.h"
#include "context.h"
#include "facility_tests.h"
#include "pyhooks.h"

#include <nlohmann/json.hpp>

#include "miso_helper.h"

namespace misoenrichment {
namespace gpr_reactor_test {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Composition::Ptr valid_composition() {
  cyclus::CompMap valid_cm;
  valid_cm[922350000] = 1.5;
  valid_cm[922380000] = 98.5;

  return cyclus::Composition::CreateFromMass(valid_cm);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Composition::Ptr invalid_composition() {
  cyclus::CompMap invalid_cm;
  invalid_cm[922350000] = 1.5;
  invalid_cm[922380000] = 98.5;
  invalid_cm[10010000] = 1;

  return cyclus::Composition::CreateFromMass(invalid_cm);
}

}  // namespace gpr_reactor_test

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Functions handling the initiation and ending of the tests.
void GprReactorTest::SetUp() {
    cyclus::PyStart();

    facility = new GprReactor(tc.get());
    InitParameters();
    SetUpGprReactor();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void GprReactorTest::InitParameters() {
    in_commods = std::vector<std::string>(1, "fresh_fuel");
    out_commods = std::vector<std::string>(1, "spent_fuel");
    in_recipes = std::vector<std::string>(1, "fresh_fuel_recipe");
    out_recipes = std::vector<std::string>(1, "spent_fuel_recipe");
    fuel_prefs = std::vector<double>(1, 1);
    n_assem_core = 2;
    n_assem_batch = 2;
    assem_size = 50000;
    n_assem_fresh = 1;
    n_assem_spent = 1;
    latitude = 0.;
    longitude = 0.;
    decom_transmute_all = false;
    cycle_time = 3;
    refuel_time = 1;
    temperature = 350;
    power_output = 2400;

    /*
    // The calculations below need to be done because of the way the SRS kernel
    // handles the power output (namely as a power density).
    const int n_assemblies_tot = 600;
    const int n_assemblies_model = 18;
    const double feet_to_m = 0.3048;  // in m * ft^-1
    double assembly_length = 12;  // in ft
    assembly_length *= feet_to_m;  // in m

    power_output *= n_assemblies_model / n_assemblies_tot / assembly_length;
    */
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void GprReactorTest::SetUpGprReactor() {
  using cyclus::Material;

  facility->in_commods = in_commods;
  facility->out_commods = out_commods;
  facility->in_recipes = in_recipes;
  facility->out_recipes = out_recipes;
  facility->fuel_prefs = fuel_prefs;
  facility->n_assem_core = n_assem_core;
  facility->n_assem_batch = n_assem_batch;
  facility->assem_size = assem_size;
  facility->n_assem_fresh = n_assem_fresh;
  facility->n_assem_spent = n_assem_spent;
  facility->latitude = latitude;
  facility->longitude = longitude;
  facility->decom_transmute_all = decom_transmute_all;
  facility->cycle_time = cycle_time;
  facility->refuel_time = refuel_time;
  facility->power_output = power_output;
  facility->temperature = temperature;

  for (int i = 0; i < n_assem_core; ++i) {
    Material::Ptr mat = Material::CreateUntracked(assem_size,
        gpr_reactor_test::valid_composition());
    facility->core.Push(mat);
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void GprReactorTest::TearDown() {
    delete facility;
    cyclus::PyStop();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Unit tests
TEST_F(GprReactorTest, ExportCompositions) {
  using cyclus::Composition;

  Composition::Ptr invalid_comp = gpr_reactor_test::invalid_composition();
  EXPECT_THROW(DoCompositionToOutFile(invalid_comp, true),
               cyclus::ValueError);

  Composition::Ptr valid_comp = gpr_reactor_test::valid_composition();
  EXPECT_NO_THROW(DoCompositionToOutFile(valid_comp, true));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(GprReactorTest, ImportCompositions) {
  using cyclus::Composition;

  Composition::Ptr valid_comp = gpr_reactor_test::valid_composition();
  cyclus::CompMap valid_cm = valid_comp->mass();

  // Prepare import composition-tests.
  nlohmann::json json_object;
  cyclus::compmath::Normalize(&valid_cm);
  for (const std::pair<int,double>& x : valid_cm) {
    // Multiply the mass fraction with the total mass of the core.
    json_object["spent_fuel_composition"][std::to_string(x.first)] =
        n_assem_core * assem_size * x.second;
  }
  std::ofstream file("gpr_reactor_spent_fuel_composition.json",
                     std::ofstream::out | std::ofstream::trunc);
  file << json_object;
  file.close();

  // Perform and test the imports
  ASSERT_NO_THROW(DoImportSpentFuelComposition(n_assem_core * assem_size));
  cyclus::CompMap return_cm = DoImportSpentFuelComposition(
      n_assem_core * assem_size)->mass();
  cyclus::compmath::Normalize(&return_cm);
  EXPECT_TRUE(cyclus::compmath::AlmostEq(valid_cm, return_cm, kEpsCompMap));

  nlohmann::json json_object2;
  json_object2["test"] = -1;
  file.open("gpr_reactor_spent_fuel_composition.json",
            std::ofstream::out | std::ofstream::trunc);
  file << json_object2;
  file.close();
  EXPECT_THROW(DoImportSpentFuelComposition(n_assem_core * assem_size / 2),
               cyclus::IOError);
  std::remove("gpr_reactor_spent_fuel_composition.json");
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(GprReactorTest, SpentFuelWaste) {
  using cyclus::Composition;

  Composition::Ptr valid_comp = gpr_reactor_test::valid_composition();
  cyclus::CompMap valid_cm = valid_comp->mass();

  // Prepare import composition-tests.
  nlohmann::json json_object;
  cyclus::compmath::Normalize(&valid_cm);
  for (const std::pair<int,double>& x : valid_cm) {
    // In this test, multiply by half a core (as opposed to the test in
    // `GprReactorTest.ImportCompositions`) because we want to test if half of
    // the spent fuel gets set to hydrogen.
    json_object["spent_fuel_composition"][std::to_string(x.first)] =
        x.second * n_assem_core * assem_size / 2.;
  }
  std::ofstream file("gpr_reactor_spent_fuel_composition.json",
                     std::ofstream::out | std::ofstream::trunc);
  file << json_object;
  file.close();

  // Perform and test the import.
  cyclus::CompMap return_cm = DoImportSpentFuelComposition(
      n_assem_core * assem_size)->mass();
  cyclus::compmath::Normalize(&return_cm);
  valid_cm[10010000] = 1.;
  cyclus::compmath::Normalize(&valid_cm);
  EXPECT_TRUE(cyclus::compmath::AlmostEq(valid_cm, return_cm, kEpsCompMap));
 
  std::remove("gpr_reactor_spent_fuel_composition.json");
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(GprReactorTest, TransmuteFuel) {
  DoTransmute();
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Below are standard tests introduced by cycstub. Probably not so useful but
// they don't do harm, either.
TEST_F(GprReactorTest, Print) {
  EXPECT_NO_THROW(std::string s = facility->str());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(GprReactorTest, Tick) {
  EXPECT_NO_THROW(facility->Tick());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(GprReactorTest, Tock) {
  EXPECT_NO_THROW(facility->Tock());
}

}  // namespace misoenrichment

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Do Not Touch! Below section required for connection with Cyclus
cyclus::Agent* GprReactorConstructor(cyclus::Context* ctx) {
  return new misoenrichment::GprReactor(ctx);
}
// Required to get functionality in cyclus agent unit tests library
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED
INSTANTIATE_TEST_CASE_P(GprReactor, FacilityTests,
                        ::testing::Values(&GprReactorConstructor));
INSTANTIATE_TEST_CASE_P(GprReactor, AgentTests,
                        ::testing::Values(&GprReactorConstructor));
