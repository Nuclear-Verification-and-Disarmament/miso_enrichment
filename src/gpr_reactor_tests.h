#ifndef MISOENRICHMENT_SRC_GPR_REACTOR_TESTS_H_
#define MISOENRICHMENT_SRC_GPR_REACTOR_TESTS_H_

#include <gtest/gtest.h>

#include "agent_tests.h"
#include "gpr_reactor.h"

namespace misoenrichment {

class GprReactorTest : public ::testing::Test {
 protected:
  cyclus::TestContext tc;
  GprReactor* facility;

  virtual void SetUp();
  virtual void TearDown();

  void InitParameters();
  void SetUpGprReactor();

  bool decom_transmute_all;

  double assem_size;
  double latitude;
  double longitude;
  double power_output;
  double temperature;

  int n_assem_core;
  int n_assem_batch;
  int n_assem_fresh;
  int n_assem_spent;
  int cycle_time;
  int refuel_time;

  std::vector<std::string> in_commods;
  std::vector<std::string> out_commods;
  std::vector<std::string> in_recipes;
  std::vector<double> fuel_prefs;

  // The Do* functions are a hack: GprReactorTest is a friend class to
  // GprReactor and it can therefore access private functions. However,
  // all of the tests are sub-classes of the fixture and they cannot access
  // the private functions, hence we need to use the Do* functions as an
  // intermediary. See, e.g., 4th bullet point here:
  // https://github.com/google/googletest/blob/master/googletest/docs/advanced.md#testing-private-code
  //
  // TODO implement functions below
  inline cyclus::Composition::Ptr DoImportSpentFuelComposition(double qty) {
    return facility->ImportSpentFuelComposition_(qty);
  }
  inline void DoCompositionToOutFile(cyclus::Composition::Ptr comp,
                                     bool delete_txt) {
    facility->CompositionToOutFile_(comp, delete_txt);
  }
  inline void DoTransmute() {
    facility->Transmute_();
  }
};

}  // namespace misoenrichment

#endif  // MISOENRICHMENT_SRC_GPR_REACTOR_TESTS_H_
