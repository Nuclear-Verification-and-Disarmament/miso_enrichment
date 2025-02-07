#ifndef MISOENRICHMENT_SRC_SOURCE_TESTS_H_
#define MISOENRICHMENT_SRC_SOURCE_TESTS_H_

#include "var_recipe_source.h"

#include <gtest/gtest.h>

#include <boost/shared_ptr.hpp>

#include "agent_tests.h"
#include "context.h"
#include "exchange_context.h"
#include "facility_tests.h"
#include "material.h"
#include "mock_sim.h"

namespace misoenrichment {

class VarRecipeSourceTest : public ::testing::Test {
 protected:
  using NucDist = std::pair<std::string, std::vector<double> >;
  using NucDistMap = std::map<int, NucDist>;
  using VarOutRecipe = std::pair<std::string, NucDistMap>;

  cyclus::TestContext tc;
  TestFacility* trader;

  VarRecipeSource* src_facility;

  double capacity;
  int simdur;
  NucDistMap nuc_map;
  VarOutRecipe var_out_recipe;

  std::string out_commod;

  boost::shared_ptr<cyclus::ExchangeContext<cyclus::Material> > GetContext(
      int nreqs, std::string commodity);

  void InitParameters();
  void SetUp();
  void SetUpVarRecipeSource();
  void TearDown();

  // Helper functions to get and set certain parameters of `VarRecipeSource`
  inline VarOutRecipe varoutrecipe(misoenrichment::VarRecipeSource* s) {
    return s->var_out_recipe_;
  }

  inline std::string outcommod(misoenrichment::VarRecipeSource* s) {
    return s->out_commod;
  }

  inline std::vector<int> throughput_times(misoenrichment::VarRecipeSource* s) {
    return s->throughput_times;
  }

  inline std::vector<double> throughput_vals(misoenrichment::VarRecipeSource* s) {
    return s->throughput_vals;
  }

  inline void update_var_out_recipe(misoenrichment::VarRecipeSource* s,
                                    VarOutRecipe v) {
    s->var_out_recipe_ = v;
  }

  inline void update_var_out_recipe(misoenrichment::VarRecipeSource* s) {
    update_var_out_recipe(s, var_out_recipe);
  }

  inline void check_var_out_recipe(misoenrichment::VarRecipeSource* s){
    s->CheckVarOutRecipe_();
  }
};

}  // namespace misoenrichment

#endif  // MISOENRICHMENT_SRC_SOURCE_TESTS_H_
