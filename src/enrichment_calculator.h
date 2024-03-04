#ifndef MISOENRICHMENT_SRC_ENRICHMENT_CALCULATOR_H_
#define MISOENRICHMENT_SRC_ENRICHMENT_CALCULATOR_H_

#include <map>
#include <vector>

#include <gtest/gtest.h>

#include "include/cppoptlib/problem.h"  // from external/CppNumericalSolvers
#include "composition.h"

class EnrichmentProblem;

namespace misoenrichment {

class EnrichmentCalculator {
 public:
  friend class EnrichmentProblem;
  FRIEND_TEST(EnrichmentCalculatorTest, AssignmentOperator);

  EnrichmentCalculator();
  EnrichmentCalculator(double gamma_235);
  EnrichmentCalculator(cyclus::Composition::Ptr feed_comp,
                       double target_product_assay,
                       double target_tails_assay, double gamma,
                       double feed_qty, double product_qty,
                       double max_swu, bool use_downblending=true,
                       bool use_integer_stages=true);
  EnrichmentCalculator(const EnrichmentCalculator& e);
  EnrichmentCalculator& operator= (const EnrichmentCalculator& e);

  void PPrint();

  void BuildMatchedAbundanceRatioCascade();

  void SetInput(cyclus::Composition::Ptr new_feed_composition,
      double new_target_product_assay, double new_target_tails_assay,
      double new_feed_qty, double new_product_qty, double new_max_swu,
      double gamma_235, bool use_downblending);

  void EnrichmentOutput(cyclus::Composition::Ptr& product_comp,
                        cyclus::Composition::Ptr& tails_comp, double& feed_used,
                        double& swu_used, double& product_produced,
                        double& tails_produced, double& n_enrich,
                        double& n_strip);
  void ProductOutput(cyclus::Composition::Ptr&, double&);

  inline double FeedUsed() { return feed_qty; }
  inline double SwuUsed() { return swu; }

 private:
  bool use_downblending;  // Use only in conjunction with `use_integer_stages`.
  bool use_integer_stages;  // Else use floating-point number of stages


  // All assays, compositions, cyclus::CompMap are assumed to be atom fraction
  cyclus::CompMap feed_composition;
  cyclus::CompMap product_composition;
  cyclus::CompMap tails_composition;

  double target_product_assay;
  double target_tails_assay;
  // Units of all of the streams are kg timestep^-1
  double target_feed_qty;
  double target_product_qty;
  double feed_qty;
  double product_qty;
  double tails_qty;
  double max_swu;  // in kg SWU timestep^-1
  double swu = 0;  // Separative work that has been performed
                   // in kg SWU timestep^-1

  const std::vector<int> isotopes;
  std::map<int,double> separation_factors;
  std::map<int,double> alpha_star;

  // Number of stages in the enriching and in the stripping section
  double n_enriching;
  double n_stripping;

  double gamma_235;  // The overall separation factor for U-235

  void CalculateGammaAlphaStar_();
  void CalculateIntegerStages_();
  void CalculateDecimalStages_();
  void CalculateFlows_();
  void CalculateSwu_();
  void CalculateConcentrations_();
  void Downblend_();
  void CalculateSums(double& sum_e, double& sum_s);

  double ValueFunction_(const cyclus::CompMap& composition);
};

class EnrichmentProblem : public cppoptlib::Problem<double> {
 public:
  EnrichmentProblem(EnrichmentCalculator*);

  // Function to be minimised (must *not* be renamed).
  // Calculates the relative difference between actual and desired uranium
  // assay for an enrichment process with a given staging; for both product and
  // tails output.
  double value(const cppoptlib::Problem<double>::TVector &staging);

 private:
  EnrichmentCalculator* calculator;
};

}  // namespace misoenrichment

#endif  // MISOENRICHMENT_SRC_ENRICHMENT_CALCULATOR_H_
