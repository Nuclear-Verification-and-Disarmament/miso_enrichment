#ifndef MISOENRICHMENT_SRC_MISO_HELPER_H_
#define MISOENRICHMENT_SRC_MISO_HELPER_H_

#include <map>
#include <vector>

#include "composition.h"
#include "material.h"

namespace misoenrichment {

namespace misotest {

bool EXPECT_TRUE_compmap(const cyclus::CompMap& cm1, 
                         const cyclus::CompMap& cm2, 
                         const double kEpsComp);

}   // namespace misotest

const double eps_compmap = 1e-10;

void IsotopesNucID(std::vector<int>& isotopes);
int IsotopeToNucID(int isotope);
int NucIDToIsotope(int nuc_id);
int ResBufIdx(
    const std::vector<cyclus::Composition::Ptr>& buf_compositions,
    const cyclus::Composition::Ptr& in_comp);

double MIsoAtomAssay(cyclus::Composition::Ptr comp);
double MIsoAtomAssay(cyclus::Material::Ptr rsrc);
double MIsoAtomAssay(std::map<int,double> compmap);

double MIsoMassAssay(cyclus::Composition::Ptr comp);
double MIsoMassAssay(cyclus::Material::Ptr rsrc);
double MIsoMassAssay(std::map<int,double> compmap);

double MIsoAtomFrac(cyclus::Composition::Ptr composition,
                            int isotope);
double MIsoAtomFrac(cyclus::Material::Ptr rsrc, int isotope);
double MIsoAtomFrac(cyclus::CompMap compmap, int isotope);

double MIsoMassFrac(cyclus::Composition::Ptr composition, 
                            int isotope);
double MIsoMassFrac(cyclus::Material::Ptr rsrc, int isotope);
double MIsoMassFrac(std::map<int,double> compmap, int isotope);

// Calculates the stage separation factor for all isotopes starting from 
// the given U235 overall separation factor.
//
// Returns a map containing the stage separation factors for all U isotopes
// with the keys being the isotopes' mass.
//
// Note that the stage separation factor is defined as the ratio of 
// abundance ratio in product to abundance ratio in tails. This method 
// follows Houston G. Wood 'Effects of Separation Processes on Minor 
// Uranium Isotopes in Enrichment Cascades'. In: Science and Global 
// Security, 16:26--36 (2008). ISSN: 0892-9882.
// DOI: 10.1080/08929880802361796
std::map<int,double> CalculateSeparationFactor(double gamma_235);

} // namespace misoenrichment

#endif  // MISOENRICHMENT_SRC_MISO_HELPER_H_
