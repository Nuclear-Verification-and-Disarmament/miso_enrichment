#ifndef MULTIISOTOPEENRICHMENT_SRC_ENRICHMENT_CALCULATOR_H_
#define MULTIISOTOPEENRICHMENT_SRC_ENRICHMENT_CALCULATOR_H_

#include <gtest/gtest.h>

#include "composition.h"
#include "multi_isotope_helper.h"

namespace multiisotopeenrichment {

class EnrichmentCalculator {
 public:
  EnrichmentCalculator() {};  //TODO this constructor can be removed later
  EnrichmentCalculator(cyclus::Composition::Ptr feed_comp, 
                       double desired_product_assay, 
                       double desired_tails_assay, double gamma);
  EnrichmentCalculator(cyclus::Composition::Ptr feed_comp,
                       double desired_product_assay,
                       double desired_tails_assay, double gamma,
                       double feed_qty=1e299, double product_qty=1e299,
                       double tails_qty=1e299, double max_swu=1e299);
  // TODO in the above constructor it might not make sense to keep the 
  // default arguments for feed_qty and for product_qty. This will be
  // determined in later steps of the implementation.

  void BuildMatchedAbundanceRatioCascade();
  void SetFeedComp(cyclus::Composition::Ptr feed_comp);
  
  void EnrichmentOutput(cyclus::CompMap& product_comp, 
                        cyclus::CompMap& tails_comp, int& n_enrich, 
                        int& n_strip);
  
  FRIEND_TEST(EnrichmentCalculatorTest, ConcentrationDifference);

 private:
  cyclus::CompMap feed_composition;
  cyclus::CompMap product_composition;
  cyclus::CompMap tails_composition;

  double design_product_assay;
  double design_tails_assay;
  
  // Units of all of the streams are kg/month
  double feed_qty;
  double product_qty;
  double tails_qty;
  double max_swu;  // in kg SWU month^-1
  double swu = 0;  // Swu that has been used, in kg SWU month^-1

  std::vector<int> isotopes;
  //IsotopesNucID(isotopes);
  std::map<int,double> separation_factors;
  std::map<int,double> alpha_star;

  // Number of stages in the enriching and in the stripping section, 
  // respectively. They are stored as double to facilitate calculations, 
  // however the values will also be whole numbers.
  double n_enriching;
  double n_stripping;
  
  double gamma_235;  // The overall separation factor for U-235
  
  void CalculateNStages(double &n_stages);
  void CalculateFlows();
  void EnrichmentOutput(
      cyclus::CompMap& product_comp, cyclus::CompMap& tails_comp, 
      double& feed_used, double& swu_used, double& product_produced, 
      double& tails_produced, int& n_enrich, int& n_strip);
  double CalculateConcentrations();
  double ConcentrationDifference();
};

}  // namespace multiisotopeenrichment

#endif  // MULTIISOTOPEENRICHMENT_SRC_ENRICHMENT_CALCULATOR_H_

