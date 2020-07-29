#ifndef MISOENRICHMENT_SRC_ENRICHMENT_CALCULATOR_H_
#define MISOENRICHMENT_SRC_ENRICHMENT_CALCULATOR_H_

#include <map>
#include <vector>

#include <gtest/gtest.h>

#include "composition.h"

namespace misoenrichment {

class EnrichmentCalculator {
 public:
  EnrichmentCalculator();
  EnrichmentCalculator(double gamma_235);
  EnrichmentCalculator(cyclus::Composition::Ptr feed_comp,
                       double target_product_assay,
                       double target_tails_assay, double gamma,
                       double feed_qty=1e299, double product_qty=1e299,
                       double max_swu=1e299);
  // TODO in the above constructor it might not make sense to keep the 
  // default arguments for feed_qty and for product_qty. This will be
  // determined in later steps of the implementation.
  EnrichmentCalculator& operator= (const EnrichmentCalculator& e);
  
  void PPrint();

  void BuildMatchedAbundanceRatioCascade();
  /*
  TODO delete new_tails_qty because it should not be the limiting factor?
  void SetInput(cyclus::Composition::Ptr new_feed_composition,
    double new_target_product_assay, double new_target_tails_assay, 
    double new_feed_qty, double new_product_qty, double new_tails_qty, 
    double new_max_swu);
  */
  void SetInput(cyclus::Composition::Ptr new_feed_composition,
      double new_target_product_assay, double new_target_tails_assay, 
      double new_feed_qty, double new_product_qty, double new_max_swu,
      double gamma_235); 

  void EnrichmentOutput(
      cyclus::CompMap& product_comp, cyclus::CompMap& tails_comp, 
      double& feed_used, double& swu_used, double& product_produced, 
      double& tails_produced, int& n_enrich, int& n_strip);
  
  inline double FeedUsed() { return feed_qty; }
  inline double SwuUsed() { return swu; }

  FRIEND_TEST(EnrichmentCalculatorTest, AssignmentOperator);

 private:

  cyclus::CompMap feed_composition;
  cyclus::CompMap product_composition;
  cyclus::CompMap tails_composition;

  double target_product_assay;
  double target_tails_assay;
  // Units of all of the streams are kg/month
  double target_feed_qty;
  double target_product_qty;
  double feed_qty;
  double product_qty;
  double tails_qty;
  double max_swu;  // in kg SWU month^-1
  double swu = 0;  // Separative work that has been performed
                   // in kg SWU month^-1

  // TODO declare vector as const?
  std::vector<int> isotopes;
  //IsotopesNucID(isotopes);
  std::map<int,double> separation_factors;
  std::map<int,double> alpha_star;

  // Number of stages in the enriching and in the stripping section, 
  // respectively. They are stored as double to facilitate calculations, 
  // however the values will also be whole numbers.
  int n_enriching;
  int n_stripping;
  
  double gamma_235;  // The overall separation factor for U-235
  
  void CalculateGammaAlphaStar_();
  void CalculateNStages_();
  void CalculateFlows_();
  void CalculateSwu_();
  void CalculateConcentrations_();
  double ValueFunction_(const cyclus::CompMap& composition);
};

}  // namespace misoenrichment

#endif  // MISOENRICHMENT_SRC_ENRICHMENT_CALCULATOR_H_

