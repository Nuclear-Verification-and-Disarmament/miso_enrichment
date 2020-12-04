#include "miso_helper.h"

#include <algorithm>
#include <iterator>
#include <sstream>

#include "comp_math.h"
#include "error.h"

namespace misoenrichment {

namespace misotest {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bool CompareCompMap(cyclus::CompMap cm1, cyclus::CompMap cm2) {
  std::vector<int> isotopes;
  IsotopesNucID(isotopes);
  std::vector<int>::iterator it;
  // The following for-loop has been added to ensure that the all of the
  // uranium keys are present in both compmaps, else the comparison fails.
  for (it = isotopes.begin(); it != isotopes.end(); it++) {
    cm1[*it] += 1e-299;
    cm2[*it] += 1e-299;
  }  

  bool result = cyclus::compmath::AlmostEq(cm1, cm2, kEpsCompMap); 
  if (!result) {
    std::cout << "Value of: cm1\n"
              << "Actual:\n";
    cyclus::CompMap::iterator it;
    for (it = cm1.begin(); it!= cm1.end(); it++) {
      std::cout << "          " << it->first << ": " << it->second << "\n";
    }
    std::cout << "Expected: \n";
    for (it = cm2.begin(); it!= cm2.end(); it++) {
      std::cout << "          " << it->first << ": " << it->second << "\n";
    }
    std::cout << "\n";
  }
  return result;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Composition::Ptr comp_depletedU() {
  cyclus::CompMap comp;
  comp[922350000] = 0.1;
  comp[922380000] = 99.9;

  return cyclus::Composition::CreateFromMass(comp);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Composition::Ptr comp_natU() {
  // TODO check where this is used and if the hydrogen could be removed
  cyclus::CompMap comp;
  comp[922340000] = 5.5e-3;
  comp[922350000] = 0.711;
  comp[922380000] = 99.2835;
  comp[10010000] = 10;  // insert hydrogen to check if it is filtered out
  
  return cyclus::Composition::CreateFromMass(comp);
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Composition::Ptr comp_weapongradeU() {
  cyclus::CompMap comp;
  comp[922340000] = 0.00780791;
  comp[922350000] = 0.91020719;
  comp[922380000] = 0.08198490;

  return cyclus::Composition::CreateFromMass(comp);
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Material::Ptr mat_natU() {
  cyclus::Composition::Ptr comp = comp_natU();
  double qty = 1.;
  return cyclus::Material::CreateUntracked(qty, comp);
};

}  // namespace misotest

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void IsotopesNucID(std::vector<int>& isotopes) {
  isotopes = {232, 233, 234, 235, 236, 238};
  for (int i = 0; i < isotopes.size(); i++) {
    isotopes[i] = (92*1000 + isotopes[i]) * 10000;
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int IsotopeToNucID(int isotope) {
  std::vector<int> isotopes = {232, 233, 234, 235, 236, 238};
  std::vector<int>::iterator it;
  
  it = std::find(isotopes.begin(), isotopes.end(), isotope);
  if (it == isotopes.end()) {
    throw cyclus::ValueError("Invalid (non-uranium) isotope!");
  }
  return (92*1000 + isotope) * 10000;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int NucIDToIsotope(int nuc_id) {
  std::vector<int> isotopes;
  IsotopesNucID(isotopes);
  std::vector<int>::iterator it;
  
  it = std::find(isotopes.begin(), isotopes.end(), nuc_id);
  if (it == isotopes.end()) {
    throw cyclus::ValueError("Invalid (non-uranium) isotope!");
  }
  return nuc_id/10000 - 92*1000;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int ResBufIdx(
    const std::vector<cyclus::Composition::Ptr>& buf_compositions,
    const cyclus::Composition::Ptr& in_comp) {
  cyclus::CompMap in_compmap = in_comp->atom();
  cyclus::compmath::Normalize(&in_compmap);

  for (const cyclus::Composition::Ptr& buf_comp : buf_compositions) {
    cyclus::CompMap buf_compmap = buf_comp->atom();
    cyclus::compmath::Normalize(&buf_compmap);
    
    if (cyclus::compmath::AlmostEq(in_compmap, buf_compmap, kEpsCompMap)) {
      int i = &buf_comp - &buf_compositions[0];
      return i;
    }
  }
  return -1;  // if element is not in buf_compositions
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoAtomAssay(cyclus::Composition::Ptr comp) {
  return MIsoAtomFrac(comp, IsotopeToNucID(235));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoAtomAssay(cyclus::Material::Ptr rsrc) {
  return MIsoAtomFrac(rsrc, IsotopeToNucID(235));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoMassAssay(cyclus::Composition::Ptr comp) {
  return MIsoMassFrac(comp, IsotopeToNucID(235));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoMassAssay(cyclus::Material::Ptr rsrc) {
  return MIsoMassFrac(rsrc, IsotopeToNucID(235));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoAtomFrac(cyclus::Composition::Ptr composition, 
                            int isotope) {
  return MIsoFrac(composition->atom(), isotope);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoAtomFrac(cyclus::Material::Ptr rsrc, int isotope) {
  return MIsoAtomFrac(rsrc->comp(), isotope);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoMassFrac(cyclus::Composition::Ptr composition, 
                            int isotope) {
  return MIsoFrac(composition->mass(), isotope);
}
  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoMassFrac(cyclus::Material::Ptr rsrc, int isotope) {
  return MIsoMassFrac(rsrc->comp(), isotope);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoAssay(cyclus::CompMap compmap) {
  return MIsoFrac(compmap, IsotopeToNucID(235));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoFrac(cyclus::CompMap compmap, int isotope) {
  std::vector<int> isotopes;
  IsotopesNucID(isotopes);
  
  double isotope_assay = 0;
  double uranium_atom_frac = 0;
  
  if (isotope < 10010000) {
    std::stringstream ss;
    ss << "Isotope id '" << isotope << "'is not a valid NucID!";
    throw cyclus::ValueError(ss.str());
  }
  // Get total uranium mole fraction, all non-uranium elements are not 
  // considered here as they are directly sent to the tails.
  for (int i : isotopes) {
    if (compmap.find(i) != compmap.end()) {
      uranium_atom_frac += compmap.at(i);
      if (i==isotope) {
        isotope_assay = compmap.at(i);
      }
    }
  }
  isotope_assay /= uranium_atom_frac;
  return isotope_assay;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::map<int,double> CalculateSeparationFactor(
    double gamma_235, std::string enrichment_method) {
  std::vector<int> isotopes;
  IsotopesNucID(isotopes);
  std::map<int,double> separation_factors;
  
  //TODO remove the cout statement below
  std::cout << "enrichment method = " << enrichment_method << "\n\n\n\n";
 
  // We consider U-238 to be the key component hence the mass differences
  // are calculated with respect to this isotope.
  if (enrichment_method == "centrifuge") {
    for (int i : isotopes) {
      // Here, the mass difference is considered hence it does not matter
      // if one considers the mass of U or of UF6.
      double delta_mass = 238. - NucIDToIsotope(i);
      separation_factors[i] = 1. + delta_mass*(gamma_235-1.) / (238.-235.);
    }
  } else if (enrichment_method == "diffusion") {
    const double hexafluoride = 6. * 19.;
    const double key_isotope_hex = 238. + hexafluoride;
    for (int i : isotopes) {
      // Here, the molecular weight of UF6 is considered.
      double mass_hex =  NucIDToIsotope(i) + hexafluoride;
      separation_factors[i] = std::pow(key_isotope_hex / mass_hex, 0.5); 
    }
  } else {
    std::string msg = "Variable 'enrichment_method' has to be either "
                      "'centrifuge' or 'diffusion'.";
    throw cyclus::ValueError(msg);
  }
  
  return separation_factors;
}

} // namespace misoenrichment
