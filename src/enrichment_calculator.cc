#include "enrichment_calculator.h"

#include <cmath>
#include <iostream>
#include <string>

#include "comp_math.h"
#include "cyc_limits.h"
#include "error.h"

#include "miso_helper.h"
#include "../CppNumericalSolvers/include/cppoptlib/solver/bfgssolver.h"

namespace misoenrichment {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Constructor delegation only possible from C++11 onwards, CMake checks if
// C++11 is supported.
EnrichmentCalculator::EnrichmentCalculator() : EnrichmentCalculator(1.4) {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Constructor delegation only possible from C++11 onwards, CMake checks if
// C++11 is supported.
EnrichmentCalculator::EnrichmentCalculator(double gamma_235) :
  EnrichmentCalculator(misotest::comp_reprocessedU(), 0.05, 0.003,
                       gamma_235, 1, 1, 1e299, true) {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
EnrichmentCalculator::EnrichmentCalculator(
    cyclus::Composition::Ptr feed_comp,
    double target_product_assay, double target_tails_assay,
    double gamma_235, double feed_qty, double product_qty,
    double max_swu, bool use_downblending, bool use_integer_stages) :
      feed_composition(feed_comp->atom()),
      target_product_assay(target_product_assay),
      target_tails_assay(target_tails_assay),
      gamma_235(gamma_235),
      target_feed_qty(feed_qty),
      target_product_qty(product_qty),
      feed_qty(0.), product_qty(0.),
      max_swu(max_swu),
      use_downblending(use_downblending),
      use_integer_stages(use_integer_stages),
      isotopes(IsotopesNucID()) {
  if (feed_qty==1e299 && product_qty==1e299 && max_swu==1e299) {
    // TODO think about whether one or two of these variables have to be
    // defined. Additionally, add an exception that should be thrown.
  }
  if (use_downblending && !use_integer_stages) {
    throw cyclus::ValueError(
      "'use_integer_stages' must be 'true' if 'use_downblending' is 'true'"
    );
  }
  cyclus::compmath::Normalize(&feed_composition);
  CalculateGammaAlphaStar_();

  BuildMatchedAbundanceRatioCascade();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
EnrichmentCalculator::EnrichmentCalculator(const EnrichmentCalculator& e) :
    feed_composition(e.feed_composition),
    product_composition(e.product_composition),
    tails_composition(e.tails_composition),
    target_product_assay(e.target_product_assay),
    target_tails_assay(e.target_tails_assay),
    target_feed_qty(e.target_feed_qty),
    target_product_qty(e.target_product_qty), max_swu(e.max_swu),
    use_downblending(e.use_downblending), isotopes(IsotopesNucID()),
    gamma_235(e.gamma_235) {
  CalculateGammaAlphaStar_();
  BuildMatchedAbundanceRatioCascade();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::CalculateGammaAlphaStar_() {
  separation_factors = CalculateSeparationFactor(gamma_235);
  for (int i : isotopes) {
    // E. von Halle Eq. (15)
    alpha_star[i] = separation_factors[i]
                    / std::sqrt(separation_factors[IsotopeToNucID(235)]);
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
EnrichmentCalculator& EnrichmentCalculator::operator= (
    const EnrichmentCalculator& e) {

  feed_composition = e.feed_composition;
  product_composition = e.product_composition;
  tails_composition = e.tails_composition;

  target_product_assay = e.target_product_assay;
  target_tails_assay = e.target_tails_assay;
  target_feed_qty = e.target_feed_qty;
  target_product_qty = e.target_product_qty;
  max_swu = e.max_swu;  // in kg SWU month^-1

  use_downblending = e.use_downblending;

  gamma_235 = e.gamma_235;
  CalculateGammaAlphaStar_();

  // TODO Check why the recalculated variables are not copied
  BuildMatchedAbundanceRatioCascade();

  return *this;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::PPrint() {
  std::cout << "- - - - - - - - - - - - - - - - - - - - - -\n"
            << "MIso Enrichment Calculator with parameters:\n"
            << "  Target product assay   " << target_product_assay << "\n"
            << "  Target tails assay     " << target_tails_assay << "\n"
            << "  Maximum SWU            " << max_swu << "\n\n"
            << "  Feed quantity          " << feed_qty << "\n"
            << "  Product quantity       " << product_qty << "\n"
            << "  Tails quantity         " << tails_qty << "\n"
            << "  Separative work used   " << swu << "\n\n"
            << "  n(enriching)           " << n_enriching << "\n"
            << "  n(stripping)           " << n_stripping << "\n"
            << "  Separation factors         232     233      234      235"
            << "      236      238\n                         ";
  for (int nuc : isotopes) {
    printf("%6.4f   ", separation_factors[nuc]);
  }
  std::cout << "\n  Compositions\n"
            << "  Isotope         Feed     Product       Tails\n";
  for (int nuc : isotopes) {
    printf("      %3d   %10.4e  %10.4e  %10.4e\n", NucIDToIsotope(nuc),
        feed_composition[nuc], product_composition[nuc],
        tails_composition[nuc]);
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::SetInput(
    cyclus::Composition::Ptr new_feed_composition,
    double new_target_product_assay, double new_target_tails_assay,
    double new_feed_qty, double new_product_qty, double new_max_swu,
    double new_gamma_235, bool new_use_downblending) {

  feed_composition = new_feed_composition->atom();
  cyclus::compmath::Normalize(&feed_composition);

  if (new_gamma_235 != gamma_235) {
    gamma_235 = new_gamma_235;
    CalculateGammaAlphaStar_();
  }
  target_product_assay = new_target_product_assay;
  target_tails_assay = new_target_tails_assay;

  target_feed_qty = new_feed_qty;
  target_product_qty = new_product_qty;
  max_swu = new_max_swu;

  use_downblending = new_use_downblending;

  // TODO Think about reimplementing this part to ensure that only the bare
  // minimum of update calculations are performed
  /*
  // If any of the concentrations change, then redesign the cascade from
  // scratch.
  if (!cyclus::compmath::AlmostEq(new_compmap, feed_composition,
                                  kEpsCompMap)
      || !cyclus::AlmostEq(new_target_product_assay, target_product_assay)
      || !cyclus::AlmostEq(new_target_tails_assay, target_tails_assay)) {
    feed_composition = new_compmap;
    target_product_assay = new_target_product_assay;
    target_tails_assay = new_target_tails_assay;

    target_feed_qty = new_feed_qty;
    target_product_qty = new_product_qty;
    max_swu = new_max_swu;

    BuildMatchedAbundanceRatioCascade();
  } else if (!cyclus::AlmostEq(new_feed_qty, target_feed_qty)
      || !cyclus::AlmostEq(new_product_qty, target_product_qty)
      || !cyclus::AlmostEq(new_max_swu, max_swu)) {
    target_feed_qty = new_feed_qty;
    target_product_qty = new_product_qty;
    max_swu = new_max_swu;

    if (feed_qty==1e299 && product_qty==1e299 && max_swu==1e299) {
      // TODO think about whether one or two of these variables have to be
      // defined. Additionally, add an exception that should be thrown.
    }
    CalculateFlows_();

    // if this code snippet is used again then add downblending!!

  }
  */
  BuildMatchedAbundanceRatioCascade();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::EnrichmentOutput(
    cyclus::Composition::Ptr& product_comp, cyclus::Composition::Ptr& tails_comp,
    double& feed_used, double& swu_used, double& product_produced,
    double& tails_produced, double& n_enrich, double& n_strip) {

  // This step is needed to prevent a 'double free' error. It is caused for
  // unknown reasons by cyclus::Material::ExtractComp and
  // cyclus::compmath::ApplyThreshold. See also Cyclus issue #1524:
  // https://github.com/cyclus/cyclus/issues/1524
  cyclus::CompMap cm;
  for (const auto& x : product_composition) {
    if (x.second > 0) {
      cm[x.first] = x.second;
    }
  }
  product_comp = cyclus::Composition::CreateFromAtom(cm);
  tails_comp = cyclus::Composition::CreateFromAtom(tails_composition);

  feed_used = feed_qty;
  swu_used = swu;
  product_produced = product_qty;
  tails_produced = tails_qty;

  n_enrich = n_enriching;
  n_strip = n_stripping;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::ProductOutput(
    cyclus::Composition::Ptr& old_product_comp, double& old_product_qty) {
  old_product_comp = cyclus::Composition::CreateFromAtom(product_composition);
  old_product_qty = product_qty;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::BuildMatchedAbundanceRatioCascade() {
  if (use_integer_stages) {
    CalculateIntegerStages_();
  } else {
    CalculateDecimalStages_();
  }
  CalculateFlows_();
  if (use_downblending) {
    Downblend_();
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::CalculateIntegerStages_() {
  // The target concentrations should always be reached or exceeded (i.e.,
  // at least equal U235 concentration in product, at most equal
  // U235 concentration in tails).
  n_enriching = 0;
  n_stripping = 0;
  do {
    n_enriching++;
    CalculateConcentrations_();
  } while (MIsoAssay(product_composition) < target_product_assay
           && n_enriching <= kIterMax);
  do {
    n_stripping++;
    CalculateConcentrations_();
  } while (MIsoAssay(tails_composition) > target_tails_assay
           && n_stripping <= kIterMax);

  if ((n_enriching == kIterMax) || (n_stripping == kIterMax)) {
    throw cyclus::Error("Unable to determine the number of stages!");
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::CalculateDecimalStages_() {
  // Solver determines the number of stages sucht that the relative difference
  // between actual and desired uranium assay in product and tails gets
  // minimised.
  EnrichmentProblem problem(this);
  cppoptlib::BfgsSolver<EnrichmentProblem> solver;

  // Offer multiple starting values to the solver to improve convergence.
  std::vector<double> n_init_stages({1, 10, 50});
  cppoptlib::Problem<double>::TVector staging(2);

  bool found_solution = false;
  for (double n_init_enriching : n_init_stages) {
    for (double n_init_stripping : n_init_stages) {
      staging[0] = n_init_enriching;
      staging[1] = n_init_stripping;
      solver.minimize(problem, staging);
      // Rough sanity checks to ensure that no non-sensical solution got used.
      // Maybe this will be replaced later by a bounded problem.
      if (
        problem(staging) < 1e-6
        && staging[0] > 0 && staging[0] < 100
        && staging[1] > 0 && staging[1] < 100
      ) {
        found_solution = true;
        break;
      }
    }
    if (found_solution) {
      break;
    }
  }
  if (!found_solution) {
    PPrint();
    throw cyclus::Error("Did not manage to determine the correct staging!");
  }
  n_enriching = staging[0];
  n_stripping = staging[1];
  CalculateConcentrations_();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::CalculateFlows_() {
  double sum_e;
  double sum_s;
  CalculateSums(sum_e, sum_s);

  // In the following, it is determined if the feed or the product quantity
  // available is a constraint. For this, the target feed and product
  // quantities are used and then it is compared which of both streams is
  // the constraining factor.
  feed_qty = target_product_qty / sum_e;  // Eq. (47)
  product_qty = target_feed_qty * sum_e;  // Eq. (47)

  // If we produce less product than desired then the feed is contraining.
  bool feed_is_constraint = product_qty < target_product_qty;
  product_qty = feed_is_constraint ? product_qty : target_product_qty;
  feed_qty = feed_is_constraint ? target_feed_qty : feed_qty;
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
double EnrichmentCalculator::ValueFunction_(
    const cyclus::CompMap& composition) {
  const int NUCID_235 = IsotopeToNucID(235);
  const int NUCID_238 = IsotopeToNucID(238);
  double value = 0.;

  try {
    composition.at(NUCID_235);
    composition.at(NUCID_238);
  } catch (const std::out_of_range& err) {
    if (composition.size()==0) {
      // This case can happen, e.g., during initalisation and is not a bug.
      return value;
    }
    // Else, throw an error.
    std::stringstream msg;
    msg << "No U-235 or U-238 present in composition passed to "
        << "'EnrichmentCalculator::ValueFunction_'.";
    throw cyclus::KeyError(msg.str());
  }

  for (int i : isotopes) {
    try {
      composition.at(i);
    } catch (const std::out_of_range& err) {
      continue;
    }
    double k = (separation_factors[i]-1)
               / (separation_factors[NUCID_235]-1);
    if (cyclus::AlmostEq(k, 0.5)) {
      // This formula is not included in  de la Garza 1963, it is taken
      // from the preceding article, see Eq. (26) in:
      // A. de la Garza et al., 'Multicomponent isotope separation in
      // cascades'. Chemical Engineering Science 15, pp. 188-209 (1961).
      value += std::log(composition.at(i) / composition.at(NUCID_238));
    } else {
      value += composition.at(i) / (2*k - 1);
    }
  }
  value *= std::log(composition.at(NUCID_235) / composition.at(NUCID_238));

  return value;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::CalculateConcentrations_() {
  // Variable naming follows E. von Halle, the equation numbers also refer
  // to his article.
  std::map<int,double> e;
  std::map<int,double> s;
  for (int i : isotopes) {
    // Eq. (37)
    e[i] = 1. / alpha_star[i]
           / (1.-std::pow(alpha_star[i], -n_enriching));
    // Eq. (39)
    s[i] = 1. / alpha_star[i]
           / (std::pow(alpha_star[i], n_stripping+1.)-1.);
  }

  double sum_e;
  double sum_s;
  CalculateSums(sum_e, sum_s);

  // Calculate the compositions of product and tails.
  for (int i : isotopes) {
    double atom_frac = MIsoFrac(feed_composition, i);
    product_composition[i] = e[i] * atom_frac / (e[i]+s[i]) / sum_e;
    tails_composition[i] = s[i] * atom_frac / (e[i]+s[i]) / sum_s;
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::Downblend_() {
  double feed_assay = MIsoAssay(feed_composition);
  double product_assay = MIsoAssay(product_composition);

  if (product_assay - target_product_assay < 0.00005) {
    return;
  }

  // Both const values store the initial value passed by MIsoEnrich while
  // the non-const target quantities are used to calculate the amount of
  // feed and product material for enrichment (excluding downblending).
  const double TARGET_FEED_QTY = target_feed_qty;
  const double TARGET_PRODUCT_QTY = target_product_qty;

  // Quantity of blending feed needed per unit of product
  double blend_feed_per_product = (product_assay-target_product_assay)
                                  / (target_product_assay-feed_assay);

  // Check whether product or feed is the constraining factor and adapt
  // the corresponding quantity.
  if (cyclus::AlmostEq(product_qty, TARGET_PRODUCT_QTY)) {
    target_product_qty /= 1 + blend_feed_per_product;
    CalculateFlows_();
  } else if (cyclus::AlmostEq(feed_qty, TARGET_FEED_QTY)) {
    double sum_e;
    double sum_s;
    CalculateSums(sum_e, sum_s);

    target_feed_qty /= 1 + blend_feed_per_product*sum_e;
    CalculateFlows_();
  }

  // Calculate the blending feed needed. If the SWU is the constraint or in
  // 'close-call' cases (e.g., product is the constraint but feed is nearly
  // completely used or vice-versa), check that both product and feed
  // limits are not exceeded.
  double blend_feed = blend_feed_per_product * product_qty;
  if (blend_feed+feed_qty > TARGET_FEED_QTY
      && !cyclus::AlmostEq(blend_feed+feed_qty, TARGET_FEED_QTY)) {
    double sum_e;
    double sum_s;
    CalculateSums(sum_e, sum_s);

    target_feed_qty /= 1 + blend_feed_per_product*sum_e;
    CalculateFlows_();

    blend_feed = blend_feed_per_product * product_qty;
  }
  if (blend_feed+product_qty > TARGET_PRODUCT_QTY
      && !cyclus::AlmostEq(blend_feed+product_qty, TARGET_PRODUCT_QTY)) {
    target_product_qty /= 1 + blend_feed_per_product;
    CalculateFlows_();

    blend_feed = blend_feed_per_product * product_qty;
  }

  // Finalise the downblending:
  // Calculate the downblended product composition.
  for (int i : isotopes) {
    product_composition[i] = (product_composition[i]*product_qty
                              + feed_composition[i]*blend_feed)
                             / (product_qty+blend_feed);
  }
  // Add blending feed to mass balances.
  feed_qty += blend_feed;
  product_qty += blend_feed;

  return;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void EnrichmentCalculator::CalculateSums(double& sum_e, double& sum_s) {
  // Variable naming follows E. von Halle, the equation numbers also refer
  // to his article.
  sum_e = 0;
  sum_s = 0;

  for (int i : isotopes) {
    double atom_frac = MIsoFrac(feed_composition, i);
    // Eq. (37)
    double e = 1. / alpha_star[i]
                  / (1-std::pow(alpha_star[i],-n_enriching));
    // Eq. (39)
    double s = 1. / alpha_star[i]
                  / (std::pow(alpha_star[i], n_stripping+1)-1);
    sum_e += e * atom_frac / (e+s);  // right-hand side of Eq. (47)
    sum_s += s * atom_frac / (e+s);  // right-hand side of Eq. (50)
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
EnrichmentProblem::EnrichmentProblem(EnrichmentCalculator* c) :
  calculator(c) {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double EnrichmentProblem::value(const cppoptlib::Problem<double>::TVector &staging) {
  calculator->n_enriching = staging[0];
  calculator->n_stripping = staging[1];
  calculator->CalculateConcentrations_();

  double target_p = calculator->target_product_assay;
  double target_t = calculator->target_tails_assay;

  return pow((calculator->product_composition[922350000] - target_p) / target_p, 2)
         + pow((calculator->tails_composition[922350000] - target_t) / target_t, 2);
}

}  // namespace misoenrichment
