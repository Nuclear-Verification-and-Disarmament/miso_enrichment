#include <cmath>
#include <map>
#include <vector>

#include "error.h"
#include "comp_math.h"

#include "enrichment_calculator.h"
#include "multi_isotope_helper.h"

namespace multiisotopeenrichment {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
EnrichmentCalculator::EnrichmentCalculator(
    cyclus::Composition::Ptr feed_comp, double desired_product_assay,
    double desired_tails_assay, double gamma) {
  feed_composition = feed_comp->atom();
  design_product_assay = desired_product_assay;
  design_tails_assay = desired_tails_assay;

  gamma_235 = gamma;
  
  BuildMatchedAbundanceRatioCascade();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::BuildMatchedAbundanceRatioCascade() {
  CalculateNStages(n_enriching);
  CalculateNStages(n_stripping);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::EnrichmentOutput(
    cyclus::CompMap& product_comp, cyclus::CompMap& tails_comp, 
    int& n_enrich, int& n_strip) {
  product_comp = product_composition;
  tails_comp = tails_composition;

  // Verify again that the numbers of stages are whole numbers.
  if (std::fmod(n_enriching, 1.) < 1e-9) {
    throw cyclus::ValueError("n_enriching is not a whole number!");
  }
  n_enrich = (int) n_enriching;
  if (std::fmod(n_stripping, 1.) < 1e-9) {
    throw cyclus::ValueError("n_stripping is not a whole number!");
  }
  n_strip = (int) n_stripping;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// TODO REMOVE THIS FUNCTION AND REPLACE IT WITH STH MORE COMPREHENSIVE
/*
void EnrichmentCalculator::SetFeedComp(cyclus::Composition::Ptr feed_comp) {
  cyclus::CompMap feed_compmap = feed_comp->atom();
  cyclus::compmath::Normalize(&feed_compmap);
  cyclus::compmath::Normalize(&feed_composition);

  // Check if feed is different from previous one, if no then do not 
  // rebuild the cascade.
  if (!cyclus::compmath::AlmostEq(feed_compmap, feed_composition, 1e-12)) {
    feed_composition = feed_comp;
    BuildMatchedAbundanceRatioCascade();
  }
}
*/


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::CalculateNStages(double &n_stages) {
  const double iter_max = 200;
  const double eps = 1e-4;
  double delta = 1e299;
  double previous_delta;
  
  n_stages = 0;
  do {
    n_stages++;
    previous_delta = delta;
    delta = CalculateConcentrations();
  } while (delta < previous_delta && n_stages != iter_max);

  if (n_stages == iter_max) {
    throw cyclus::Error("Unable to determine the number of stages!");
  }
  
  n_stages = n_stages - (1+eps);
  n_stages = previous_delta > CalculateConcentrations()
             ? n_stages + eps : n_stages + (1+eps);

  // ensure that the number of stages is an integer (stored as double)
  if (std::fmod(n_stages, 1.) < 1e-9) {
    throw cyclus::ValueError("n_stages is not a whole number!");
  }
  CalculateConcentrations();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double EnrichmentCalculator::CalculateConcentrations() {
  // Calculates gamma = alpha*beta for all isotopes.
  // Variable naming follows E. von Halle, the equation numbers also refer
  // to his article.
  std::vector<int> isotopes;
  IsotopesNucID(isotopes);
  std::map<int,double> separation_factors;
  std::map<int,double> alpha_star;
  std::map<int,double> e;
  std::map<int,double> s;
  
  double alpha_235;
  double atom_frac;
  double delta_concentration;
  double e_sum = 0;
  double s_sum = 0;
  
  // Define the above declared variables.
  // TODO: Calculate alpha_235 of centrifuge such that alpha = beta
  separation_factors = CalculateSeparationFactor(alpha_235);
  for (int i : isotopes) {
    // Eq. (15)
    alpha_star[i] = separation_factors[i]
                    / std::sqrt(separation_factors[IsotopeToNucID(235)]); 
    // Eq. (37)
    e[i] = 1. / alpha_star[i] 
           / (1-std::pow(alpha_star[i],-n_enriching));
    // Eq. (39)
    s[i] = 1. / alpha_star[i] 
           / (std::pow(alpha_star[i], n_stripping+1)-1);

    atom_frac = MultiIsotopeAtomFrac(feed_composition, i);
    e_sum += e[i] * atom_frac / (e[i]+s[i]);  // Eq. (48) denominator
    s_sum += s[i] * atom_frac / (e[i]+s[i]);  // Eq. (51) denominator
  }
  
  // Calculate the compositions of product and tails.
  for (int i : isotopes) {
    atom_frac = MultiIsotopeAtomFrac(feed_composition, i);
    product_composition[i] = e[i] * atom_frac / (e[i]+s[i]) / e_sum;
    tails_composition[i] = s[i] * atom_frac / (e[i]+s[i]) / s_sum;
  }

  delta_concentration = ConcentrationDifference();
  return delta_concentration;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double EnrichmentCalculator::ConcentrationDifference() {
  int nuc_id235 = IsotopeToNucID(235);
  double product_assay = product_composition[nuc_id235];
  double tails_assay = tails_composition[nuc_id235];

  double delta_product = (product_assay-design_product_assay)
                         / design_product_assay;
  double delta_tails = (tails_assay-design_tails_assay) 
                       / design_tails_assay;
  double delta = std::sqrt(std::pow(delta_product, 2.) 
                           + std::pow(delta_tails, 2.));
  return delta;
}

}  // namespace multiisotopeenrichment
