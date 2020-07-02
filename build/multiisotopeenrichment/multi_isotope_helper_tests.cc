
#include <gtest/gtest.h>

#include "composition.h"
#include "cyc_limits.h"
#include "material.h"
#include "pyne.h"

#include "multi_isotope_helper.h"

namespace multiisotopeenrichment {

namespace multiisotopehelpertest {

cyclus::Composition::Ptr comp_natU() {
  std::map<int,double> m;
  m[922340000] = 5.5e-3;
  m[922350000] = 0.711;
  m[922380000] = 99.2835;
  m[10010000] = 10;  // insert hydrogen to check if it is filtered out
  return cyclus::Composition::CreateFromMass(m);
};
cyclus::Material::Ptr mat_natU() {
  cyclus::Composition::Ptr comp = comp_natU();
  double qty = 1.;
  return cyclus::Material::CreateUntracked(qty, comp);
};

} // namespace multiisotopehelpertest

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST(MultiIsotopeHelperTest, NucIDConversion) {
  std::vector<int> isotopes;
  IsotopesNucID(isotopes);
  
  for (int i : isotopes) {
    int isotope = NucIDToIsotope(i);
    EXPECT_EQ(IsotopeToNucID(isotope), i);
  }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST(MultiIsotopeHelperTest, CheckFractionsComposition) {
  double expected_mass235 = 0.00711;
  double expected_atom235 = 0.00711 / pyne::atomic_mass(922350000)
                            / (5.5e-5/pyne::atomic_mass(922340000)
                               + 0.00711/pyne::atomic_mass(922350000)
                               + 0.992835/pyne::atomic_mass(922380000));

  cyclus::Composition::Ptr comp = multiisotopehelpertest::comp_natU();
  EXPECT_TRUE(cyclus::AlmostEq(MultiIsotopeAtomAssay(comp), 
                               expected_atom235));
  EXPECT_TRUE(cyclus::AlmostEq(MultiIsotopeAtomFrac(comp, 922350000),
                               expected_atom235));
  EXPECT_TRUE(cyclus::AlmostEq(MultiIsotopeMassAssay(comp), 
                               expected_mass235));
  EXPECT_TRUE(cyclus::AlmostEq(MultiIsotopeMassFrac(comp, 922350000), 
                               expected_mass235));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST(MultiIsotopeHelperTest, CheckFractionsMaterial) {
  double expected_mass235 = 0.00711;
  double expected_atom235 = 0.00711 / pyne::atomic_mass(922350000)
                            / (5.5e-5/pyne::atomic_mass(922340000)
                               + 0.00711/pyne::atomic_mass(922350000)
                               + 0.992835/pyne::atomic_mass(922380000));

  cyclus::Material::Ptr mat = multiisotopehelpertest::mat_natU();
  EXPECT_TRUE(cyclus::AlmostEq(MultiIsotopeAtomAssay(mat), 
                               expected_atom235));
  EXPECT_TRUE(cyclus::AlmostEq(MultiIsotopeAtomFrac(mat, 922350000),
                               expected_atom235));
  EXPECT_TRUE(cyclus::AlmostEq(MultiIsotopeMassAssay(mat), 
                               expected_mass235));
  EXPECT_TRUE(cyclus::AlmostEq(MultiIsotopeMassFrac(mat, 922350000), 
                               expected_mass235));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST(MultiIsotopeHelperTest, SeparationFactor) {
  // The values below are valid for alpha*beta = 1.6
  double alpha_235 = std::pow(1.6, 0.5);
  std::map<int,double> separation_factor = CalculateSeparationFactor(
                                                              alpha_235);

  std::vector<int> isotopes;
  IsotopesNucID(isotopes);
  std::map<int,double> expected;
  expected[922320000] = 2.2;
  expected[922330000] = 2.0;
  expected[922340000] = 1.8;
  expected[922350000] = 1.6;
  expected[922360000] = 1.4;
  expected[922380000] = 1.0;
  
  for (int i : isotopes) {
    EXPECT_TRUE(cyclus::AlmostEq(separation_factor[i], expected[i]));
  }
}

} // namespace multiisotopeenrichment

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED

/*
cyclus::Agent* EnrichmentConstructor(cyclus::Context* ctx) {
  return new cycamore::Enrichment(ctx);
}

// required to get functionality in cyclus agent unit tests library
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INSTANTIATE_TEST_CASE_P(EnrichmentFac, FacilityTests,
                        Values(&EnrichmentConstructor));
INSTANTIATE_TEST_CASE_P(EnrichmentFac, AgentTests,
                        Values(&EnrichmentConstructor));
*/
