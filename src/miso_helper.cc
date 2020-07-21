#include "miso_helper.h"

#include <algorithm>
#include <iterator>

#include "comp_math.h"
#include "error.h"

namespace misoenrichment {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void IsotopesNucID(std::vector<int> &isotopes) {
  isotopes = {232, 233, 234, 235, 236, 238};
  for (int i = 0; i < isotopes.size(); i++) {
    isotopes[i] = (92*1000 + isotopes[i]) * 10000;
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int IsotopeToNucID(int isotope) {
  std::vector<int> isotopes = {232, 233, 234, 235, 236, 238};
  std::vector<int>::iterator it;
  
  it = std::find(isotopes.begin(), isotopes.end(), isotope);
  if (it == isotopes.end()) {
    throw cyclus::ValueError("Invalid (non-uranium) isotope!");
  }
  return (92*1000 + isotope) * 10000;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int ResBufIdx(
    const std::vector<cyclus::Composition::Ptr>& buf_compositions,
    const cyclus::Composition::Ptr& in_comp) {
  cyclus::CompMap in_compmap = in_comp->atom();
  cyclus::compmath::Normalize(&in_compmap);

  for (const cyclus::Composition::Ptr& buf_comp : buf_compositions) {
    cyclus::CompMap buf_compmap = buf_comp->atom();
    cyclus::compmath::Normalize(&buf_compmap);
    
    if (cyclus::compmath::AlmostEq(in_compmap, buf_compmap, eps_compmap)) {
      int i = &buf_comp - &buf_compositions[0];
      return i;
    }
  }
  return buf_compositions.size();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoAtomAssay(cyclus::Composition::Ptr comp) {
  return MIsoAtomFrac(comp, IsotopeToNucID(235));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoAtomAssay(cyclus::Material::Ptr rsrc) {
  return MIsoAtomFrac(rsrc, IsotopeToNucID(235));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoAtomAssay(std::map<int,double> compmap) {
  return MIsoAtomFrac(compmap, IsotopeToNucID(235));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoMassAssay(cyclus::Composition::Ptr comp) {
  return MIsoMassFrac(comp, IsotopeToNucID(235));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoMassAssay(cyclus::Material::Ptr rsrc) {
  return MIsoMassFrac(rsrc, IsotopeToNucID(235));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoMassAssay(std::map<int,double> compmap) {
  return MIsoMassFrac(compmap, IsotopeToNucID(235));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoAtomFrac(cyclus::Composition::Ptr composition, 
                            int isotope) {
  return MIsoAtomFrac(composition->atom(), isotope);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoAtomFrac(cyclus::Material::Ptr rsrc, int isotope) {
  return MIsoAtomFrac(rsrc->comp(), isotope);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoAtomFrac(cyclus::CompMap compmap, int isotope) {
  std::vector<int> isotopes;
  IsotopesNucID(isotopes);
  
  double isotope_assay;
  double uranium_atom_frac = 0;

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

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoMassFrac(cyclus::Composition::Ptr composition, 
                            int isotope) {
  return MIsoMassFrac(composition->mass(), isotope);
}
  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoMassFrac(cyclus::Material::Ptr rsrc, int isotope) {
  return MIsoMassFrac(rsrc->comp(), isotope);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double MIsoMassFrac(std::map<int,double> compmap, int isotope) {
  std::vector<int> isotopes;
  IsotopesNucID(isotopes);

  double isotope_assay;
  double uranium_mass_frac = 0;

  // Get total uranium mass fraction, all non-uranium elements are not 
  // considered here as they are directly sent to the tails.
  for (int i : isotopes) {
    if (compmap.find(i) != compmap.end()) {
      uranium_mass_frac += compmap.at(i);
      if (i==isotope) {
        isotope_assay = compmap.at(i);
      }
    }
  }
  isotope_assay /= uranium_mass_frac;
  return isotope_assay;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::map<int,double> CalculateSeparationFactor(double gamma_235) {
  std::vector<int> isotopes;
  IsotopesNucID(isotopes);
  std::map<int,double> separation_factors;
  
  // We consider U-238 to be the key component hence the mass differences
  // are calculated with respect to this isotope.
  for (int i : isotopes) {
    double delta_mass = 238. - NucIDToIsotope(i);
    double gamma = 1. + delta_mass*(gamma_235-1.) / (238.-235.);
    separation_factors[i] = gamma;
  }
  return separation_factors;
}

} // namespace misoenrichment
