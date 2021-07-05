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

cyclus::Composition::Ptr valid_composition() {
  cyclus::CompMap valid_cm;
  valid_cm[922350000] = 1.5;
  valid_cm[922380000] = 98.5;

  return cyclus::Composition::CreateFromMass(valid_cm);
}

cyclus::Composition::Ptr invalid_composition() {
  cyclus::CompMap invalid_cm;
  invalid_cm[922350000] = 1.5;
  invalid_cm[922380000] = 98.5;
  invalid_cm[10010000] = 1;

  return cyclus::Composition::CreateFromMass(invalid_cm);
}

// The compositions below are taken from CNERG's cycamore module, see
// https://github.com/cyclus/cycamore/blob/master/src/reactor_tests.cc
// Copyright under BSD-3 License, University of Wisconsin Computational Nuclear
// Engineering Research Group, 2010-2016.
cyclus::Composition::Ptr c_uox() {
  cyclus::CompMap m;
  m[pyne::nucname::id("u235")] = 0.04;
  m[pyne::nucname::id("u238")] = 0.96;
  return cyclus::Composition::CreateFromMass(m);
};

cyclus::Composition::Ptr c_mox() {
  cyclus::CompMap m;
  m[pyne::nucname::id("u235")] = .7;
  m[pyne::nucname::id("u238")] = 100;
  m[pyne::nucname::id("pu239")] = 3.3;
  return cyclus::Composition::CreateFromMass(m);
};

cyclus::Composition::Ptr c_spentuox() {
  cyclus::CompMap m;
  m[pyne::nucname::id("u235")] =  .8;
  m[pyne::nucname::id("u238")] =  100;
  m[pyne::nucname::id("pu239")] = 1;
  return cyclus::Composition::CreateFromMass(m);
};

cyclus::Composition::Ptr c_spentmox() {
  cyclus::CompMap m;
  m[pyne::nucname::id("u235")] =  .2;
  m[pyne::nucname::id("u238")] =  100;
  m[pyne::nucname::id("pu239")] = .9;
  return cyclus::Composition::CreateFromMass(m);
};

cyclus::Composition::Ptr c_water() {
  cyclus::CompMap m;
  m[pyne::nucname::id("O16")] =  1;
  m[pyne::nucname::id("H1")] =  2;
  return cyclus::Composition::CreateFromAtom(m);
};

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
    fuel_prefs = std::vector<double>(1, 1);
    n_assem_core = 1;
    n_assem_batch = 1;
    assem_size = 110820;
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
    cyclus::Material::Ptr mat = cyclus::Material::CreateUntracked(assem_size,
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
    // Convert the nuc id to a name (e.g., 922350001 --> U235M) to match the
    // output of the `spentfuelgpr` program and multiply the mass fraction with
    // the total mass of the core.
    std::string nuc_name = pyne::nucname::name(x.first);
    double mass = n_assem_core * assem_size * x.second;
    json_object["spent_fuel_composition"][nuc_name] = mass;
  }
  std::ofstream file(DoGetInFname(), std::ofstream::out | std::ofstream::trunc);
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
  file.open(DoGetInFname(), std::ofstream::out | std::ofstream::trunc);
  file << json_object2;
  file.close();
  EXPECT_THROW(DoImportSpentFuelComposition(n_assem_core * assem_size / 2),
               cyclus::IOError);
  std::remove(DoGetInFname().c_str());
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
    // Convert the nuc id to a name (e.g., 922350001 --> U235M) to match the
    // output of the `spentfuelgpr` program and multiply the mass fraction with
    // the total mass of the core.
    // In this test, multiply by half a core (as opposed to the test in
    // `GprReactorTest.ImportCompositions`) because we want to test if half of
    // the spent fuel gets set to hydrogen.
    std::string nuc_name = pyne::nucname::name(x.first);
    double mass = n_assem_core * assem_size * x.second / 2.;
    json_object["spent_fuel_composition"][nuc_name] = mass;
  }
  std::ofstream file(DoGetInFname(), std::ofstream::out | std::ofstream::trunc);
  file << json_object;
  file.close();

  // Perform and test the import.
  cyclus::CompMap return_cm = DoImportSpentFuelComposition(
      n_assem_core * assem_size)->mass();
  cyclus::compmath::Normalize(&return_cm);
  valid_cm[10010000] = 1.;
  cyclus::compmath::Normalize(&valid_cm);
  EXPECT_TRUE(cyclus::compmath::AlmostEq(valid_cm, return_cm, kEpsCompMap));
 
  std::remove(DoGetInFname().c_str());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(GprReactorTest, TransmuteFuel) {
  DoTransmute();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Below are unit tests taken from CNERG's cycamore module, see
// https://github.com/cyclus/cycamore/blob/master/src/reactor_tests.cc
// Copyright under BSD-3 License, University of Wisconsin Computational Nuclear
// Engineering Research Group, 2010-2016.

// Test that with a zero refuel_time and a zero capacity fresh fuel buffer
// (the default), fuel can be ordered and the cycle started with no time step
// delay.
TEST_F(GprReactorTest, JustInTimeOrdering) {
  std::string config =
     "  <fuel_inrecipes>  <val>lwr_fresh</val>  </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>lwr_spent</val>  </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>enriched_u</val> </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>      </fuel_outcommods>  "
     "  <fuel_prefs>      <val>1.0</val>        </fuel_prefs>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <assem_size>300</assem_size>  "
     "  <n_assem_core>1</n_assem_core>  "
     "  <n_assem_batch>1</n_assem_batch>  ";

  int simdur = 50;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  sim.AddSource("enriched_u").Finalize();
  sim.AddRecipe("lwr_fresh", gpr_reactor_test::c_uox());
  sim.AddRecipe("lwr_spent", gpr_reactor_test::c_spentuox());
  int id = sim.Run();

  cyclus::QueryResult qr = sim.db().Query("Transactions", NULL);
  EXPECT_EQ(simdur, qr.rows.size()) << "failed to order+run on fresh fuel "
                                       "inside 1 time step";
}

// tests that the correct number of assemblies are popped from the core each
// cycle.
TEST_F(GprReactorTest, BatchSizes) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <assem_size>1</assem_size>  "
     "  <n_assem_core>7</n_assem_core>  "
     "  <n_assem_batch>3</n_assem_batch>  ";

  int simdur = 50;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", gpr_reactor_test::c_uox());
  sim.AddRecipe("spentuox", gpr_reactor_test::c_spentuox());
  int id = sim.Run();

  cyclus::QueryResult qr = sim.db().Query("Transactions", NULL);
  // 7 for initial core, 3 per time step for each new batch for remainder
  EXPECT_EQ(7+3*(simdur-1), qr.rows.size());
}

// tests that the refueling period between cycle end and start of the next
// cycle is honored.
TEST_F(GprReactorTest, RefuelTimes) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>4</cycle_time>  "
     "  <refuel_time>3</refuel_time>  "
     "  <assem_size>1</assem_size>  "
     "  <n_assem_core>1</n_assem_core>  "
     "  <n_assem_batch>1</n_assem_batch>  ";

  int simdur = 49;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", gpr_reactor_test::c_uox());
  sim.AddRecipe("spentuox", gpr_reactor_test::c_spentuox());
  int id = sim.Run();

  cyclus::QueryResult qr = sim.db().Query("Transactions", NULL);
  int cyclet = 4;
  int refuelt = 3;
  int n_assem_want = simdur/(cyclet+refuelt)+1; // +1 for initial core
  EXPECT_EQ(n_assem_want, qr.rows.size());
}


// tests that a reactor decommissions on time without producing
// power at the end of its lifetime.
TEST_F(GprReactorTest, DecomTimes) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>2</cycle_time>  "
     "  <refuel_time>2</refuel_time>  "
     "  <assem_size>1</assem_size>  "
     "  <n_assem_core>3</n_assem_core>  "
     "  <power_cap>1000</power_cap>  "
     "  <n_assem_batch>1</n_assem_batch>  ";

  int simdur = 12;
  int lifetime = 7;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur,
                      lifetime);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", gpr_reactor_test::c_uox());
  sim.AddRecipe("spentuox", gpr_reactor_test::c_spentuox());
  int id = sim.Run();

  // operating for 2+2 months and shutdown for 2+1
  int on_time = 4;
  std::vector<cyclus::Cond> conds;
  conds.push_back(cyclus::Cond("Value", "==", 1000));
  cyclus::QueryResult qr = sim.db().Query("TimeSeriesPower", &conds);
  EXPECT_EQ(on_time, qr.rows.size());

  int off_time = 3;
  conds.clear();
  conds.push_back(cyclus::Cond("Value", "==", 0));
  qr = sim.db().Query("TimeSeriesPower", &conds);
  EXPECT_EQ(off_time, qr.rows.size());
}


// Tests if a reactor produces power at the time of its decommission
// given a refuel_time of zero.
TEST_F(GprReactorTest, DecomZeroRefuel) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>2</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <assem_size>1</assem_size>  "
     "  <n_assem_core>3</n_assem_core>  "
     "  <power_cap>1000</power_cap>  "
     "  <n_assem_batch>1</n_assem_batch>  ";

  int simdur = 8;
  int lifetime = 6;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur,
                      lifetime);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", gpr_reactor_test::c_uox());
  sim.AddRecipe("spentuox", gpr_reactor_test::c_spentuox());
  int id = sim.Run();

  // operating for 2+2 months and shutdown for 2+1
  int on_time = 6;
  std::vector<cyclus::Cond> conds;
  conds.push_back(cyclus::Cond("Value", "==", 1000));
  cyclus::QueryResult qr = sim.db().Query("TimeSeriesPower", &conds);
  EXPECT_EQ(on_time, qr.rows.size());
}

// tests that new fuel is ordered immediately following cycle end - at the
// start of the refueling period - not before and not after. - thie is subtly
// different than RefuelTimes test and is not a duplicate of it.
TEST_F(GprReactorTest, OrderAtRefuelStart) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>4</cycle_time>  "
     "  <refuel_time>3</refuel_time>  "
     "  <assem_size>1</assem_size>  "
     "  <n_assem_core>1</n_assem_core>  "
     "  <n_assem_batch>1</n_assem_batch>  ";

  int simdur = 7;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", gpr_reactor_test::c_uox());
  sim.AddRecipe("spentuox", gpr_reactor_test::c_spentuox());
  int id = sim.Run();

  cyclus::QueryResult qr = sim.db().Query("Transactions", NULL);
  int cyclet = 4;
  int refuelt = 3;
  int n_assem_want = simdur/(cyclet+refuelt)+1; // +1 for initial core
  EXPECT_EQ(n_assem_want, qr.rows.size());
}

// tests that the reactor handles requesting multiple types of fuel correctly
// - with proper inventory constraint honoring, etc.
TEST_F(GprReactorTest, MultiFuelMix) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      <val>mox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> <val>spentmox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      <val>mox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <assem_size>1</assem_size>  "
     "  <n_assem_fresh>3</n_assem_fresh>  "
     "  <n_assem_core>3</n_assem_core>  "
     "  <n_assem_batch>3</n_assem_batch>  ";

  // it is important that the sources have cumulative capacity greater than
  // the reactor can take on a single time step - to test that inventory
  // capacity constraints are being set properly.  It is also important that
  // each source is smaller capacity thatn the reactor orders on each time
  // step to make it easy to compute+check the number of transactions.
  int simdur = 50;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  sim.AddSource("uox").capacity(2).Finalize();
  sim.AddSource("mox").capacity(2).Finalize();
  sim.AddRecipe("uox", gpr_reactor_test::c_uox());
  sim.AddRecipe("spentuox", gpr_reactor_test::c_spentuox());
  sim.AddRecipe("mox", gpr_reactor_test::c_spentuox());
  sim.AddRecipe("spentmox", gpr_reactor_test::c_spentuox());
  int id = sim.Run();

  cyclus::QueryResult qr = sim.db().Query("Transactions", NULL);
  // +3 is for fresh fuel inventory
  EXPECT_EQ(3*simdur+3, qr.rows.size());
}

// tests that the reactor halts operation when it has no more room in its
// spent fuel inventory buffer.
TEST_F(GprReactorTest, FullSpentInventory) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <assem_size>1</assem_size>  "
     "  <n_assem_core>1</n_assem_core>  "
     "  <n_assem_batch>1</n_assem_batch>  "
     "  <n_assem_spent>3</n_assem_spent>  ";

  int simdur = 10;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", gpr_reactor_test::c_uox());
  sim.AddRecipe("spentuox", gpr_reactor_test::c_spentuox());
  int id = sim.Run();

  cyclus::QueryResult qr = sim.db().Query("Transactions", NULL);
  int n_assem_spent = 3;

  // +1 is for the assembly in the core + the three in spent
  EXPECT_EQ(n_assem_spent+1, qr.rows.size());
}

// tests that the reactor shuts down, ie., does not generate power, when the
// spent fuel inventory is full and the core cannot be unloaded.
TEST_F(GprReactorTest, FullSpentInventoryShutdown) {
  std::string config =
    " <fuel_inrecipes> <val>uox</val> </fuel_inrecipes> "
    " <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes> "
    " <fuel_incommods> <val>uox</val> </fuel_incommods> "
    " <fuel_outcommods> <val>waste</val> </fuel_outcommods> "
    ""
    " <cycle_time>1</cycle_time> "
    " <refuel_time>0</refuel_time> "
    " <assem_size>1</assem_size> "
    " <n_assem_core>1</n_assem_core> "
    " <n_assem_batch>1</n_assem_batch> "
    " <n_assem_spent>1</n_assem_spent> "
    " <power_cap>100</power_cap> ";

  int simdur = 3;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", gpr_reactor_test::c_uox());
  sim.AddRecipe("spentuox", gpr_reactor_test::c_spentuox());
  int id = sim.Run();

  cyclus::QueryResult qr = sim.db().Query("TimeSeriesPower", NULL);
  EXPECT_EQ(0, qr.GetVal<double>("Value", simdur - 1));

}

// tests that the reactor cycle is delayed as expected when it is unable to
// acquire fuel in time for the next cycle start.  This checks that after a
// cycle is delayed past an original scheduled start time, as soon as enough
// fuel is received, a new cycle pattern is established starting from the
// delayed start time.
TEST_F(GprReactorTest, FuelShortage) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>7</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <assem_size>1</assem_size>  "
     "  <n_assem_core>3</n_assem_core>  "
     "  <n_assem_batch>3</n_assem_batch>  ";

  int simdur = 50;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  // provide initial full batch
  sim.AddSource("uox").lifetime(1).Finalize();
  // provide partial batch post cycle-end
  sim.AddSource("uox").start(9).lifetime(1).capacity(2).Finalize();
  // provide remainder of batch much later
  sim.AddSource("uox").start(15).Finalize();
  sim.AddRecipe("uox", gpr_reactor_test::c_uox());
  sim.AddRecipe("spentuox", gpr_reactor_test::c_spentuox());
  int id = sim.Run();

  // check that we never got a full refueled batch during refuel period
  std::vector<cyclus::Cond> conds;
  conds.push_back(cyclus::Cond("Time", "<", 15));
  cyclus::QueryResult qr = sim.db().Query("Transactions", &conds);
  EXPECT_EQ(5, qr.rows.size());

  // after being delayed past original scheduled start of new cycle, we got
  // final assembly for core.
  conds.clear();
  conds.push_back(cyclus::Cond("Time", "==", 15));
  qr = sim.db().Query("Transactions", &conds);
  EXPECT_EQ(1, qr.rows.size());

  // all during the next (delayed) cycle we shouldn't be requesting any new fuel
  conds.clear();
  conds.push_back(cyclus::Cond("Time", "<", 21));
  qr = sim.db().Query("Transactions", &conds);
  EXPECT_EQ(6, qr.rows.size());

  // as soon as this delayed cycle ends, we should be requesting/getting 3 new
  // batches
  conds.clear();
  conds.push_back(cyclus::Cond("Time", "==", 22));
  qr = sim.db().Query("Transactions", &conds);
  EXPECT_EQ(3, qr.rows.size());
}

// tests that discharged fuel is transmuted properly immediately at cycle end.
TEST_F(GprReactorTest, DischargedFuelTransmute) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>4</cycle_time>  "
     "  <refuel_time>3</refuel_time>  "
     "  <assem_size>1</assem_size>  "
     "  <n_assem_core>1</n_assem_core>  "
     "  <n_assem_batch>1</n_assem_batch>  ";

  int simdur = 7;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  sim.AddSource("uox").Finalize();
  sim.AddSink("waste").Finalize();
  sim.AddRecipe("uox", gpr_reactor_test::c_uox());
  Composition::Ptr spentuox = gpr_reactor_test::c_spentuox();
  sim.AddRecipe("spentuox", spentuox);
  int id = sim.Run();

  std::vector<cyclus::Cond> conds;
  conds.push_back(cyclus::Cond("SenderId", "==", id));
  int resid = sim.db().Query("Transactions", &conds).GetVal<int>("ResourceId");
  cyclus::Material::Ptr m = sim.GetMaterial(resid);
  cyclus::toolkit::MatQuery mq(m);
  EXPECT_EQ(spentuox->id(), m->comp()->id());
  EXPECT_TRUE(mq.mass(942390000) > 0) << "transmuted spent fuel doesn't have Pu239";
}

// tests that spent fuel is offerred on correct commods according to the
// incommod it was received on - esp when dealing with multiple fuel commods
// simultaneously.
TEST_F(GprReactorTest, SpentFuelProperCommodTracking) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      <val>mox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> <val>spentmox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      <val>mox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste1</val>   <val>waste2</val>   </fuel_outcommods>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <assem_size>1</assem_size>  "
     "  <n_assem_core>3</n_assem_core>  "
     "  <n_assem_batch>3</n_assem_batch>  ";

  int simdur = 7;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  sim.AddSource("uox").capacity(1).Finalize();
  sim.AddSource("mox").capacity(2).Finalize();
  sim.AddSink("waste1").Finalize();
  sim.AddSink("waste2").Finalize();
  sim.AddRecipe("uox", gpr_reactor_test::c_uox());
  sim.AddRecipe("spentuox", gpr_reactor_test::c_spentuox());
  sim.AddRecipe("mox", gpr_reactor_test::c_mox());
  sim.AddRecipe("spentmox", gpr_reactor_test::c_spentmox());
  int id = sim.Run();

  std::vector<cyclus::Cond> conds;
  conds.push_back(cyclus::Cond("SenderId", "==", id));
  conds.push_back(cyclus::Cond("Commodity", "==", std::string("waste1")));
  cyclus::QueryResult qr = sim.db().Query("Transactions", &conds);
  EXPECT_EQ(simdur-1, qr.rows.size());

  conds[1] = cyclus::Cond("Commodity", "==", std::string("waste2"));
  qr = sim.db().Query("Transactions", &conds);
  EXPECT_EQ(2*(simdur-1), qr.rows.size());
}

// The user can optionally omit fuel preferences.  In the case where
// preferences are adjusted, the ommitted preference vector must be populated
// with default values - if it wasn't then preferences won't be adjusted
// correctly and the reactor could segfault.  Check that this doesn't happen.
TEST_F(GprReactorTest, PrefChange) {
  // it is important that the fuel_prefs not be present in the config below.
  std::string config =
     "  <fuel_inrecipes>  <val>lwr_fresh</val>  </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>lwr_spent</val>  </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>enriched_u</val> </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>      </fuel_outcommods>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <assem_size>300</assem_size>  "
     "  <n_assem_core>1</n_assem_core>  "
     "  <n_assem_batch>1</n_assem_batch>  "
     ""
     "  <pref_change_times>   <val>25</val>         </pref_change_times>"
     "  <pref_change_commods> <val>enriched_u</val> </pref_change_commods>"
     "  <pref_change_values>  <val>-1</val>         </pref_change_values>";

  int simdur = 50;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  sim.AddSource("enriched_u").Finalize();
  sim.AddRecipe("lwr_fresh", gpr_reactor_test::c_uox());
  sim.AddRecipe("lwr_spent", gpr_reactor_test::c_spentuox());
  int id = sim.Run();

  cyclus::QueryResult qr = sim.db().Query("Transactions", NULL);
  EXPECT_EQ(25, qr.rows.size()) << "failed to adjust preferences properly";
}

TEST_F(GprReactorTest, RecipeChange) {
  // it is important that the fuel_prefs not be present in the config below.
  std::string config =
     "  <fuel_inrecipes>  <val>lwr_fresh</val>  </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>lwr_spent</val>  </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>enriched_u</val> </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>      </fuel_outcommods>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <assem_size>300</assem_size>  "
     "  <n_assem_core>1</n_assem_core>  "
     "  <n_assem_batch>1</n_assem_batch>  "
     ""
     "  <recipe_change_times>   <val>25</val>         <val>35</val>         </recipe_change_times>"
     "  <recipe_change_commods> <val>enriched_u</val> <val>enriched_u</val> </recipe_change_commods>"
     "  <recipe_change_in>      <val>water</val>      <val>water</val>      </recipe_change_in>"
     "  <recipe_change_out>     <val>lwr_spent</val>  <val>water</val>      </recipe_change_out>";

  int simdur = 50;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  sim.AddSource("enriched_u").Finalize();
  sim.AddSink("waste").Finalize();
  sim.AddRecipe("lwr_fresh", gpr_reactor_test::c_uox());
  sim.AddRecipe("lwr_spent", gpr_reactor_test::c_spentuox());
  sim.AddRecipe("water", gpr_reactor_test::c_water());
  int aid = sim.Run();

  std::vector<cyclus::Cond> conds;
  cyclus::QueryResult qr;

  // check that received recipe is not water
  conds.clear();
  conds.push_back(cyclus::Cond("Time", "==", 24));
  conds.push_back(cyclus::Cond("ReceiverId", "==", aid));
  qr = sim.db().Query("Transactions", &conds);
  cyclus::toolkit::MatQuery mq = cyclus::toolkit::MatQuery(
      sim.GetMaterial(qr.GetVal<int>("ResourceId")));

  EXPECT_TRUE(0 < mq.qty());
  EXPECT_TRUE(0 == mq.mass(pyne::nucname::id("H1")));

  // check that received recipe changed to water
  conds.clear();
  conds.push_back(cyclus::Cond("Time", "==", 26));
  conds.push_back(cyclus::Cond("ReceiverId", "==", aid));
  qr = sim.db().Query("Transactions", &conds);
  mq = cyclus::toolkit::MatQuery(sim.GetMaterial(qr.GetVal<int>("ResourceId")));

  EXPECT_TRUE(0 < mq.qty());
  EXPECT_TRUE(0 < mq.mass(pyne::nucname::id("H1")));

  // check that sent recipe is not water
  conds.clear();
  conds.push_back(cyclus::Cond("Time", "==", 34));
  conds.push_back(cyclus::Cond("SenderId", "==", aid));
  qr = sim.db().Query("Transactions", &conds);
  mq = cyclus::toolkit::MatQuery(sim.GetMaterial(qr.GetVal<int>("ResourceId")));

  EXPECT_TRUE(0 < mq.qty());
  EXPECT_TRUE(0 == mq.mass(pyne::nucname::id("H1")));

  // check that sent recipe changed to water
  conds.clear();
  conds.push_back(cyclus::Cond("Time", "==", 36));
  conds.push_back(cyclus::Cond("SenderId", "==", aid));
  qr = sim.db().Query("Transactions", &conds);
  mq = cyclus::toolkit::MatQuery(sim.GetMaterial(qr.GetVal<int>("ResourceId")));

  EXPECT_TRUE(0 < mq.qty());
  EXPECT_TRUE(0 < mq.mass(pyne::nucname::id("H1")));
}

TEST_F(GprReactorTest, Retire) {
  std::string config =
     "  <fuel_inrecipes>  <val>lwr_fresh</val>  </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>lwr_spent</val>  </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>enriched_u</val> </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>      </fuel_outcommods>  "
     ""
     "  <cycle_time>7</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <assem_size>300</assem_size>  "
     "  <n_assem_fresh>1</n_assem_fresh>  "
     "  <n_assem_core>3</n_assem_core>  "
     "  <n_assem_batch>1</n_assem_batch>  "
     "  <power_cap>1</power_cap>  "
     "";

  int dur = 50;
  int life = 36;
  int cycle_time = 7;
  int refuel_time = 0;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, dur,
                      life);
  sim.AddSource("enriched_u").Finalize();
  sim.AddSink("waste").Finalize();
  sim.AddRecipe("lwr_fresh", gpr_reactor_test::c_uox());
  sim.AddRecipe("lwr_spent", gpr_reactor_test::c_spentuox());
  int id = sim.Run();

  int ncore = 3;
  int nbatch = 1;

  // reactor should stop requesting new fresh fuel as it approaches retirement
  int nassem_recv =
      static_cast<int>(ceil(static_cast<double>(life) / 7.0)) * nbatch +
      (ncore - nbatch);

  std::vector<cyclus::Cond> conds;
  conds.push_back(cyclus::Cond("ReceiverId", "==", id));
  cyclus::QueryResult qr = sim.db().Query("Transactions", &conds);
  EXPECT_EQ(nassem_recv, qr.rows.size())
      << "failed to stop ordering near retirement";

  // reactor should discharge all fuel before/by retirement
  conds.clear();
  conds.push_back(cyclus::Cond("SenderId", "==", id));
  qr = sim.db().Query("Transactions", &conds);
  EXPECT_EQ(nassem_recv, qr.rows.size())
      << "failed to discharge all material by retirement time";

  // reactor should record power entry on the time step it retires if operating
  int time_online = life / (cycle_time + refuel_time) * cycle_time
                    + std::min(life % (cycle_time + refuel_time), cycle_time);
  conds.clear();
  conds.push_back(cyclus::Cond("AgentId", "==", id));
  conds.push_back(cyclus::Cond("Value", ">", 0));
  qr = sim.db().Query("TimeSeriesPower", &conds);
  EXPECT_EQ(time_online, qr.rows.size())
      << "failed to generate power for the correct number of time steps";
}

TEST_F(GprReactorTest, PositionInitialize) {
  std::string config =
     "  <fuel_inrecipes>  <val>lwr_fresh</val>  </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>lwr_spent</val>  </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>enriched_u</val> </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>      </fuel_outcommods>  "
     "  <fuel_prefs>      <val>1.0</val>        </fuel_prefs>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <assem_size>300</assem_size>  "
     "  <n_assem_core>1</n_assem_core>  "
     "  <n_assem_batch>1</n_assem_batch>  ";

  int simdur = 50;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  sim.AddSource("enriched_u").Finalize();
  sim.AddRecipe("lwr_fresh", gpr_reactor_test::c_uox());
  sim.AddRecipe("lwr_spent", gpr_reactor_test::c_spentuox());
  int id = sim.Run();

  cyclus::QueryResult qr = sim.db().Query("AgentPosition", NULL);
  EXPECT_EQ(qr.GetVal<double>("Latitude"), 0.0);
  EXPECT_EQ(qr.GetVal<double>("Longitude"), 0.0);
}

TEST_F(GprReactorTest, PositionInitialize2) {
  std::string config =
     "  <fuel_inrecipes>  <val>lwr_fresh</val>  </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>lwr_spent</val>  </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>enriched_u</val> </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>      </fuel_outcommods>  "
     "  <fuel_prefs>      <val>1.0</val>        </fuel_prefs>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <assem_size>300</assem_size>  "
     "  <n_assem_core>1</n_assem_core>  "
     "  <n_assem_batch>1</n_assem_batch>  "
     "  <longitude>30.0</longitude>  "
     "  <latitude>30.0</latitude>  ";

  int simdur = 50;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  sim.AddSource("enriched_u").Finalize();
  sim.AddRecipe("lwr_fresh", gpr_reactor_test::c_uox());
  sim.AddRecipe("lwr_spent", gpr_reactor_test::c_spentuox());
  int id = sim.Run();

  cyclus::QueryResult qr = sim.db().Query("AgentPosition", NULL);
  EXPECT_EQ(qr.GetVal<double>("Latitude"), 30.0);
  EXPECT_EQ(qr.GetVal<double>("Longitude"), 30.0);
}

TEST_F(GprReactorTest, ByProduct) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>1</refuel_time>  "
     "  <assem_size>1</assem_size>  "
     "  <n_assem_core>7</n_assem_core>  "
     "  <n_assem_batch>3</n_assem_batch>  "
     ""
     "  <side_products> <val>process_heat</val> </side_products>"
     "  <side_product_quantity> <val>10</val> </side_product_quantity>";

  int simdur = 10;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", gpr_reactor_test::c_uox());
  sim.AddRecipe("spentuox", gpr_reactor_test::c_spentuox());
  int id = sim.Run();

  std::vector<cyclus::Cond> conds;
  // test if it produces side products only when reactor is running
  int quantity = 10;
  conds.push_back(cyclus::Cond("Value", "==", quantity));
  cyclus::QueryResult qr = sim.db().Query("ReactorSideProducts", &conds);
  EXPECT_EQ(5, qr.rows.size());

  // test if it doesn't produce side products when reactor is refueling
  conds.clear();
    conds.push_back(cyclus::Cond("Value", "==", 0));
  qr = sim.db().Query("ReactorSideProducts", &conds);
  EXPECT_EQ(5, qr.rows.size());
}

TEST_F(GprReactorTest, MultipleByProduct) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>1</refuel_time>  "
     "  <assem_size>1</assem_size>  "
     "  <n_assem_core>7</n_assem_core>  "
     "  <n_assem_batch>3</n_assem_batch>  "
     ""
     "  <side_products> <val>process_heat</val> <val>water</val> </side_products>"
     "  <side_product_quantity> <val>10</val> <val>100</val> </side_product_quantity>";

  int simdur = 10;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:Reactor"), config, simdur);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", gpr_reactor_test::c_uox());
  sim.AddRecipe("spentuox", gpr_reactor_test::c_spentuox());
  int id = sim.Run();


  std::vector<cyclus::Cond> conds;
  // test if it produces heat when reactor is running
  int quantity = 10;
  conds.push_back(cyclus::Cond("Product", "==", std::string("process_heat")));
  conds.push_back(cyclus::Cond("Value", "==", quantity));
  cyclus::QueryResult qr = sim.db().Query("ReactorSideProducts", &conds);
  EXPECT_EQ(5, qr.rows.size());

  // test if it produces water when reactor is running
  conds.clear();
  quantity = 100;
  conds.push_back(cyclus::Cond("Product", "==", std::string("water")));
  conds.push_back(cyclus::Cond("Value", "==", quantity));
  qr = sim.db().Query("ReactorSideProducts", &conds);
  EXPECT_EQ(5, qr.rows.size());

  conds.clear();
  conds.push_back(cyclus::Cond("Value", "==", 0));
  qr = sim.db().Query("ReactorSideProducts", &conds);
  EXPECT_EQ(10, qr.rows.size());
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
