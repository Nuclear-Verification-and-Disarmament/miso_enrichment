#include <gtest/gtest.h>

#include "composition.h"
#include "material.h"
#include "pyne.h"

#include "miso_helper.h"

namespace misoenrichment {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST(MIsoHelperTest, ChooseCorrectResBuf) {
  using cyclus::Composition;
  std::vector<Composition::Ptr> comp_vec;
  comp_vec.push_back(misotest::comp_natU());
  comp_vec.push_back(misotest::comp_weapongradeU());
  
  cyclus::CompMap plutonium_cm;
  plutonium_cm[942390000] = 1.;
  Composition::Ptr plutonium = Composition::CreateFromAtom(plutonium_cm);

  EXPECT_EQ(ResBufIdx(comp_vec, misotest::comp_natU()), 0);
  EXPECT_EQ(ResBufIdx(comp_vec, plutonium), -1);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST(MIsoHelperTest, NucIDConversion) {
  std::vector<int> isotopes;
  IsotopesNucID(isotopes);
  
  for (int i : isotopes) {
    int isotope = NucIDToIsotope(i);
    EXPECT_EQ(IsotopeToNucID(isotope), i);
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST(MIsoHelperTest, CheckFractionsComposition) {
  double expected_mass235 = 0.00711;
  double expected_atom235 = 0.00711 / pyne::atomic_mass(922350000)
                            / (5.5e-5/pyne::atomic_mass(922340000)
                               + 0.00711/pyne::atomic_mass(922350000)
                               + 0.992835/pyne::atomic_mass(922380000));

  cyclus::Composition::Ptr comp = misotest::comp_natU();
  EXPECT_DOUBLE_EQ(MIsoAtomAssay(comp), expected_atom235);
  EXPECT_DOUBLE_EQ(MIsoAtomFrac(comp, 922350000), expected_atom235);
  EXPECT_DOUBLE_EQ(MIsoMassAssay(comp), expected_mass235);
  EXPECT_DOUBLE_EQ(MIsoMassFrac(comp, 922350000), expected_mass235);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST(MIsoHelperTest, CheckFractionsMaterial) {
  double expected_mass235 = 0.00711;
  double expected_atom235 = 0.00711 / pyne::atomic_mass(922350000)
                            / (5.5e-5/pyne::atomic_mass(922340000)
                               + 0.00711/pyne::atomic_mass(922350000)
                               + 0.992835/pyne::atomic_mass(922380000));

  cyclus::Material::Ptr mat = misotest::mat_natU();
  EXPECT_DOUBLE_EQ(MIsoAtomAssay(mat), expected_atom235);
  EXPECT_DOUBLE_EQ(MIsoAtomFrac(mat, 922350000), expected_atom235);
  EXPECT_DOUBLE_EQ(MIsoMassAssay(mat), expected_mass235);
  EXPECT_DOUBLE_EQ(MIsoMassFrac(mat, 922350000), expected_mass235);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST(MIsoHelperTest, SeparationFactor) {
  double gamma_235 = 1.6;
  std::string enrichment_methods[2] = {"centrifuge", "diffusion"};
  std::map<std::string,std::map<int,double> > expected;
  
  expected["centrifuge"][922320000] = 2.2;
  expected["centrifuge"][922330000] = 2.0;
  expected["centrifuge"][922340000] = 1.8;
  expected["centrifuge"][922350000] = 1.6;
  expected["centrifuge"][922360000] = 1.4;
  expected["centrifuge"][922380000] = 1.0;
  
  expected["diffusion"][922320000] = 1.008633;
  expected["diffusion"][922330000] = 1.007179;
  expected["diffusion"][922340000] = 1.005731;
  expected["diffusion"][922350000] = 1.004289;
  expected["diffusion"][922360000] = 1.002853;
  expected["diffusion"][922380000] = 1.0;
  
  std::vector<int> isotopes;
  IsotopesNucID(isotopes); 
  for (std::string method : enrichment_methods) {
    std::map<int,double> separation_factor = CalculateSeparationFactor(
        gamma_235, method);
    for (int i : isotopes) {
      // Use EXPECT_NEAR because of the limited precision of the values in
      // the diffusion case as indicated in Houston Wood's article.
      EXPECT_NEAR(separation_factor[i], expected[method][i], 1e-6);
    }
  }
}

} // namespace misoenrichment

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// required to get functionality in cyclus agent unit tests library
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED
