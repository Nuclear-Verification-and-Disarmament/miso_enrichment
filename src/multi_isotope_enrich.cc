#include "multi_isotope_enrich.h"

#include <cmath>
#include <sstream> 

namespace multiisotopeenrichment {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MultiIsotopeEnrich::MultiIsotopeEnrich(cyclus::Context* ctx)
    : cyclus::Facility(ctx),
      tails_assay(0),
      swu_capacity(0),
      max_enrich(1),
      initial_feed(0),
      feed_commod(""),
      feed_recipe(""),
      product_commod(""),
      tails_commod(""),
      order_prefs(true),
      latitude(0.0),
      longitude(0.0),
      coordinates(latitude, longitude) {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MultiIsotopeEnrich::~MultiIsotopeEnrich() {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::string MultiIsotopeEnrich::str() {
  std::stringstream ss;
  ss << cyclus::Facility::str() << " with enrichment facility parameters:";
  // TODO complete stringstream

  return ss.str();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Enrichment::Build(cyclus::Agent* parent) {
  Facility::Build(cyclus::Agent* parent);
  if (initial_feed > 0) {
    cyclus::Composition::Ptr initial_feed_comp = context()->GetRecipe(
      feed_recipe);
    int inventory_idx = ResBufIdx(feed_inv_comp, initial_feed_comp);
    feed_inv[i].Push(cyclus::Material::Create(this, initial_feed, 
                                              initial_feed_comp));
  }

  LOG(cyclus::LEV_DEBUG2, "MIsoEn") << "Multi-Isotope Enrichment Facility "
                                    << "entering the simulation: ";
  LOG(cyclus::LEV_DEBUG2, "MIsoEn") << str();
  RecordPosition();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MultiIsotopeEnrich::Tick() {
  current_swu_capacity = max_swu_capacity;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MultiIsotopeEnrich::Tock() {
  using cyclus::toolkit::RecordTimeSeries;
  
  LOG(cyclus::LEV_INFO4, "MIsoEn") << prototype() << " used "
                                   << intra_timestep_swu << " SWU"; 
  RecordTimeSeries<cyclus::toolkit::ENRICH_SWU>(this, intra_timestep_swu);
  
  LOG(cyclus::LEV_INFO4, "EnrFac") << prototype() << " used "
                                   << intra_timestep_feed << " feed";
  RecordTimeSeries<cyclus::toolkit::ENRICH_FEED>(this, 
                                                 intra_timestep_feed);
  RecordTimeSeries<double>("demand"+feed_commod, this, 
                           intra_timestep_feed);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr> 
MultiIsotopeEnrich::GetMatlRequests() {
  using cyclus::Material;
  using cyclus::RequestPortfolio;
  using cyclus::Request;

  std::set<RequestPortfolio<Material>::Ptr> ports;
  RequestPortfolio<Material>::Ptr port(new RequestPortfolio<Material>());
  Material::Ptr mat = Request();

  if (mat->quantity() > cyclus::eps_rsrc()) {
    //TODO use multiple feed commodities?
    port->AddRequest(mat, this, feed_commod);
    ports.insert(port);
  }
  return ports;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Material::Ptr Enrichment::Request() {
  double qty = std::max(0.0, feed_inv[current_feed_inv].capacity()
                             - feed_inv[current_feed_inv].quantity());
  return cyclus::Material::CreateUntracked(qty,
                                           feed_inv_comp[current_feed_inv]);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr> Enrichment::GetMatlBdis(
    cyclus::CommodMap<cyclus::Material>::type& out_requests) {
  using cyclus::Bid;
  using cyclus::BidPortfolio;
  using cyclus::Material;
  using cyclus::toolkit::RecordTimeSeries
  std::set<BidPortfolio<Material>::Ptr> ports;
  
  RecordTimeSeries<double>("supply" + tails_commod, this, 
                           tails.quantity());
                           
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// WARNING! Do not change the following this function!!! This enables your
// archetype to be dynamically loaded and any alterations will cause your
// archetype to fail.
extern "C" cyclus::Agent* ConstructMultiIsotopeEnrich(cyclus::Context* ctx) {
  return new MultiIsotopeEnrich(ctx);
}

}  // namespace multiisotopeenrichment
