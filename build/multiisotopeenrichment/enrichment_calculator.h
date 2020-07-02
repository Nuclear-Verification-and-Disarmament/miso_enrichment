
#line 1 "/Users/test/Uni/Masterarbeit/multi_isotope_enrichment/src/enrichment_calculator.h"
#ifndef MULTIISOTOPEENRICHMENT_SRC_ENRICHMENT_CALCULATOR_H_
#define MULTIISOTOPEENRICHMENT_SRC_ENRICHMENT_CALCULATOR_H_

#include <gtest/gtest.h>

#include "composition.h"

namespace multiisotopeenrichment {

class EnrichmentCalculator {
 public:
  EnrichmentCalculator() {};  //TODO this constructor can be removed later
  EnrichmentCalculator(cyclus::Composition::Ptr feed_comp, 
                       double desired_product_assay, 
                       double desired_tails_assay, double gamma);
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
  
  // Number of stages in the enriching and in the stripping section, 
  // respectively. They are stored as double to facilitate calculations, 
  // however the values will also be whole numbers.
  double n_enriching;
  double n_stripping;
  
  double gamma_235;  // The overall separation factor for U-235
  
  void CalculateNStages(double &n_stages);
  double CalculateConcentrations();
  double ConcentrationDifference();
};

}  // namespace multiisotopeenrichment

#endif  // MULTIISOTOPEENRICHMENT_SRC_ENRICHMENT_CALCULATOR_H_