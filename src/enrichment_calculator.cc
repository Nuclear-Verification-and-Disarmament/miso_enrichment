#include "enrichment_calculator.h"

#include <cmath>
#include <string>

#include "comp_math.h"
#include "cyc_limits.h"
#include "error.h"

#include "multi_isotope_helper.h"

namespace multiisotopeenrichment {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
EnrichmentCalculator::EnrichmentCalculator(
    cyclus::Composition::Ptr feed_composition, 
    double target_product_assay, double target_tails_assay, 
    double gamma_235, double feed_qty, double product_qty, 
    double max_swu) : 
      feed_composition(feed_composition->atom()),
      target_product_assay(target_product_assay),
      target_tails_assay(target_tails_assay),
      gamma_235(gamma_235), 
      feed_qty(feed_qty),
      product_qty(product_qty),
      max_swu(max_swu) {
  if (feed_qty==1e299 && product_qty==1e299 && max_swu==1e299) {
    // TODO think about whether one or two of these variables have to be 
    // defined. Additionally, add an exception that should be thrown.
  }
  cyclus::compmath::Normalize(&this->feed_composition);

  IsotopesNucID(isotopes);
  separation_factors = CalculateSeparationFactor(gamma_235);
  for (int i : isotopes) {
    // E. von Halle Eq. (15)
    alpha_star[i] = separation_factors[i]
                    / std::sqrt(separation_factors[IsotopeToNucID(235)]); 
  }
   
  BuildMatchedAbundanceRatioCascade();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
EnrichmentCalculator& EnrichmentCalculator::operator= (
    const EnrichmentCalculator& e) {
  
  feed_composition = e.feed_composition;
  product_composition = e.product_composition;
  tails_composition = e.tails_composition;

  target_product_assay = e.target_product_assay;
  target_tails_assay = e.target_tails_assay;
  
  feed_qty = e.feed_qty;
  product_qty = e.product_qty;
  tails_qty = e.tails_qty;
  max_swu = e.max_swu;  // in kg SWU month^-1
  swu = e.swu;
  
  IsotopesNucID(isotopes);
  gamma_235 = e.gamma_235;
  separation_factors = CalculateSeparationFactor(gamma_235);
  for (int i : isotopes) {
    // E. von Halle Eq. (15)
    alpha_star[i] = separation_factors[i]
                    / std::sqrt(separation_factors[IsotopeToNucID(235)]); 
  }
  
  // TODO Check why the recalculated variables are not copied
  BuildMatchedAbundanceRatioCascade();

  return *this;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::SetInput(
    cyclus::Composition::Ptr new_feed_composition,
    double new_target_product_assay, double new_target_tails_assay, 
    double new_feed_qty, double new_product_qty, double new_max_swu) {
  
  // This temporary variable is needed because feed_comp->atom() returns
  // a const cyclus::CompMap object and hence it cannot be normalised.
  // However, the normalisation is needed to be able to compare it to the
  // current feed_composition
  cyclus::CompMap new_compmap = new_feed_composition->atom();
  cyclus::compmath::Normalize(&new_compmap);

  // If any of the concentrations change, then redesign the cascade from
  // scratch.
  if (!cyclus::compmath::AlmostEq(new_compmap, feed_composition, kEpsComp)
      || !cyclus::AlmostEq(new_target_product_assay, target_product_assay)
      || !cyclus::AlmostEq(new_target_tails_assay, target_tails_assay)) {
    feed_composition = new_compmap;
    target_product_assay = new_target_product_assay;
    target_tails_assay = new_target_tails_assay;

    BuildMatchedAbundanceRatioCascade();
  }
  // If any of the flows change, then recalculate the flows. The
  // concentrations remain unaffected of this change.
  if (!cyclus::AlmostEq(new_feed_qty, feed_qty) 
      || !cyclus::AlmostEq(new_product_qty, product_qty)
      || !cyclus::AlmostEq(new_max_swu, max_swu)) {
    feed_qty = new_feed_qty;
    product_qty = new_product_qty;
    
    if (feed_qty==1e299 && product_qty==1e299 && max_swu==1e299) {
      // TODO think about whether one or two of these variables have to be
      // defined. Additionally, add an exception that should be thrown.
    }
    CalculateFlows_();
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::EnrichmentOutput(
    cyclus::CompMap& product_comp, cyclus::CompMap& tails_comp,
    double& feed_used, double& swu_used, double& product_produced, 
    double& tails_produced, int& n_enrich, int& n_strip) {

  product_comp = product_composition;
  tails_comp = tails_composition;
  
  feed_used = feed_qty;
  swu_used = swu;
  product_produced = product_qty;
  tails_produced = tails_qty;

  // Verify again that the numbers of stages are whole numbers.
  if (std::fmod(n_enriching, 1.) > 1e-9) {
    throw cyclus::ValueError("n_enriching is not a whole number!");
  }
  n_enrich = (int) n_enriching;
  if (std::fmod(n_stripping, 1.) > 1e-9) {
    throw cyclus::ValueError("n_stripping is not a whole number!");
  }
  n_strip = (int) n_stripping;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::BuildMatchedAbundanceRatioCascade() {
  CalculateNStages_(n_enriching);
  CalculateNStages_(n_stripping);

  CalculateFlows_();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::CalculateNStages_(double &n_stages) {
  const double kIterMax = 200;
  const double kEps = 1e-4;
  double delta = 1e299;
  double previous_delta;
  
  n_stages = 0;
  do {
    n_stages++;
    previous_delta = delta;
    delta = CalculateConcentrations_();
  } while (delta < previous_delta && n_stages != kIterMax);
  
  if (n_stages == kIterMax) {
    throw cyclus::Error("Unable to determine the number of stages!");
  }

  n_stages = n_stages - (1+kEps);
  n_stages = previous_delta > CalculateConcentrations_()
             ? n_stages + kEps : n_stages + (1+kEps);

  // ensure that the number of stages is an integer (stored as double)
  if (std::fmod(n_stages, 1.) > 1e-9) {
    throw cyclus::ValueError("n_stages is not a whole number!");
  }
  CalculateConcentrations_();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::CalculateFlows_() {
  double real_feed;
  double real_product;
  double real_swu;
  double sum_e = 0;
  double sum_s = 0;

  for (int i : isotopes) {
    double e;
    double s;  
    double atom_frac = MultiIsotopeAtomFrac(feed_composition, i);

    // Eq. (37)
    e = 1. / alpha_star[i] / (1-std::pow(alpha_star[i],-n_enriching));
    // Eq. (39)
    s = 1. / alpha_star[i] / (std::pow(alpha_star[i], n_stripping+1)-1);
    sum_e += e * atom_frac / (e+s);  // right-hand side of Eq. (47)
    sum_s += s * atom_frac / (e+s);  // right-hand side of Eq. (50)
  }

  // In the following, it is determined if the feed or the product quantity
  // available is a constraint. The feed used to produce `product_qty` of 
  // product is calculatedas well as the product produced by using 
  // `feed_qty` of feed. In order to differentiate, the names of these 
  // hypothetical quantities are preceded by `real_`.
  real_feed = product_qty / sum_e;  // Eq. (47)
  real_product = feed_qty * sum_e;  // Eq. (47)
  
  // If we produce less product than desired then the feed is contraining.
  bool feed_is_constraint = real_product < product_qty;
  product_qty = feed_is_constraint ? real_product : product_qty;
  feed_qty = feed_is_constraint ? feed_qty : real_feed;
  tails_qty = feed_qty * sum_s;  // Eq. (50)
  
  // Having determined the enrichment flows, calculate the separative work
  // that is performed and check if it does not exceed the SWU capacity.
  // If it does, recalculate using the given SWU capacity.
  CalculateSwu_();
  if (swu > max_swu) {
    swu = max_swu;
    
    feed_qty = swu / (ValueFunction_(product_composition)*sum_e
                      + ValueFunction_(tails_composition)*sum_s
                      - ValueFunction_(feed_composition));
    product_qty = feed_qty * sum_e;
    tails_qty = feed_qty * sum_s;
  } 
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::CalculateSwu_() {
  double v_f = ValueFunction_(feed_composition);
  double v_p = ValueFunction_(product_composition);
  double v_t = ValueFunction_(tails_composition);

  swu = v_p*product_qty + v_t*tails_qty - v_f*feed_qty;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//TODO define ValueFunction
double EnrichmentCalculator::ValueFunction_(
    const cyclus::CompMap& composition) {
  double value = 0;
  
  std::string msg("MIsoEn: Value Function not yet implemented\n");
  cyclus::Warn<cyclus::Warnings::VALIDATION_WARNING>(msg);

  /*
  std::vector<int>::iterator it;
  for (it = isotopes.begin(); it != isotopes.end(); it++) {
    
  }
  */
  return value;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double EnrichmentCalculator::CalculateConcentrations_() {
  // Variable naming follows E. von Halle, the equation numbers also refer
  // to his article.
  std::map<int,double> e;
  std::map<int,double> s;
  
  double alpha_235;
  double atom_frac;
  double delta_concentration;
  double e_sum = 0;
  double s_sum = 0;
  
  // Define the above declared variables.
  for (int i : isotopes) {
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

  delta_concentration = ConcentrationDifference_();
  return delta_concentration;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double EnrichmentCalculator::ConcentrationDifference_() {
  int nuc_id235 = IsotopeToNucID(235);
  double product_assay = product_composition[nuc_id235];
  double tails_assay = tails_composition[nuc_id235];

  double delta_product = (product_assay-target_product_assay)
                         / target_product_assay;
  double delta_tails = (tails_assay-target_tails_assay) 
                       / target_tails_assay;
  double delta = std::sqrt(std::pow(delta_product, 2.) 
                           + std::pow(delta_tails, 2.));
  return delta;
}

}  // namespace multiisotopeenrichment
