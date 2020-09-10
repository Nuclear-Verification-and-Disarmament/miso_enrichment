#include "miso_enrich.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>
#include <sstream> 
#include <vector>

#include "toolkit/timeseries.h"

#include "miso_helper.h"

namespace misoenrichment {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MIsoEnrich::MIsoEnrich(cyclus::Context* ctx)
    : cyclus::Facility(ctx),
      tails_assay(0.003),
      max_enrich(0.95),
      initial_feed(0),
      max_feed_inventory(1e299),
      feed_commod(""),
      feed_recipe(""),
      product_commod(""),
      tails_commod(""),
      order_prefs(true),
      swu_capacity(1e299),
      intra_timestep_swu(0),
      intra_timestep_feed(0),
      latitude(0.0),
      longitude(0.0),
      coordinates(latitude, longitude),
      use_downblending(true) {
  
  enrichment_calc = EnrichmentCalculator(gamma_235);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MIsoEnrich::~MIsoEnrich() {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::string MIsoEnrich::str() {
  std::stringstream ss;
  ss << cyclus::Facility::str() << " with enrichment facility parameters:";
  // TODO complete stringstream

  return ss.str();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrich::EnterNotify() {
  using cyclus::Material;
  
  cyclus::Facility::EnterNotify();

  if (swu_capacity_times[0]==-1) {
    swu_flexible = FlexibleInput<double>(this, swu_capacity_vals);
  } else {
    swu_flexible = FlexibleInput<double>(this, swu_capacity_vals, 
                                         swu_capacity_times);
  }

  if (initial_feed > 0) {
    Material::Ptr mat = Material::Create(
        this, initial_feed, context()->GetRecipe(feed_recipe));
    AddFeedMat_(mat);
  } else {
    feed_inv.push_back(cyclus::toolkit::ResBuf<cyclus::Material>());
    feed_inv.back().capacity(max_feed_inventory);
    feed_inv_comp.push_back(context()->GetRecipe(feed_recipe));
    feed_idx = 0;  // set current feed idx to the only existing inventory
  }

  LOG(cyclus::LEV_DEBUG2, "MIsoEn") << "Multi-Isotope Enrichment Facility "
                                    << "entering the simulation: ";
  LOG(cyclus::LEV_DEBUG2, "MIsoEn") << str();
  RecordPosition();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrich::AddFeedMat_(cyclus::Material::Ptr mat) {
  cyclus::Composition::Ptr comp = mat->comp();
  int push_idx = ResBufIdx(feed_inv_comp, comp);
  
  // Either directly try pushing material to the right feed inventory or 
  // create a corresponding feed inventory and add it to the vector.
  if (push_idx != -1) {
    LOG(cyclus::LEV_INFO5, "MIsoEn") << prototype() 
                                     << " is initially holding "
                                     << feed_inv[push_idx].quantity() 
                                     << " of feed in inventory no. "
                                     << push_idx << ".";
    try {  
      feed_inv[push_idx].Push(mat);
    } catch (cyclus::Error& e) {
      e.msg(Agent::InformErrorMsg(e.msg()));
    throw e;
    }
    LOG(cyclus::LEV_INFO5, "MIsoEn") << prototype() << " added " 
                                     << mat->quantity() << " of " 
                                     << feed_commod 
                                     << " to its inventory no. " << push_idx
                                     << " which is now holding " 
                                     << feed_inv[push_idx].quantity();
  } else {
    LOG(cyclus::LEV_INFO5, "MIsoEn") << prototype() << " is initially"
                                     << " holding no feed of this"
                                     << " composition and creates a new"
                                     << " inventory.";
  
    feed_inv.push_back(cyclus::toolkit::ResBuf<cyclus::Material>());
    feed_inv.back().capacity(max_feed_inventory);
    try {
      feed_inv.back().Push(mat);
    } catch (cyclus::Error& e) {
      e.msg(Agent::InformErrorMsg(e.msg()));
    }
    feed_inv_comp.push_back(comp);
    // '-1' because of index starting at 0
    feed_idx = std::distance(feed_inv.begin(), feed_inv.end()) - 1;

    LOG(cyclus::LEV_INFO5, "MIsoEn") << prototype() << " added " 
                                     << mat->quantity() << " of " 
                                     << feed_commod 
                                     << " to its new inventory (no. "
                                     << feed_idx << ") which is "
                                     << "now holding " 
                                     << feed_inv[feed_idx].quantity();  
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrich::Tick() {
  // For an unknown reason, 'UpdateValue' has to be called with a copy of
  // the 'this' pointer. When directly usiong 'this', the address passed to
  // the function is increased by 8 bits resulting later on in a
  // segmentation fault.
  cyclus::Agent* copy_ptr;
  cyclus::Agent* source_ptr = this;
  std::memcpy((void*) &copy_ptr, (void*) &source_ptr, sizeof(cyclus::Agent*));

  swu_capacity = swu_flexible.UpdateValue(copy_ptr);
  current_swu_capacity = swu_capacity;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrich::Tock() {
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
MIsoEnrich::GetMatlRequests() {
  using cyclus::Material;
  using cyclus::RequestPortfolio;
  using cyclus::Request;
  
  std::set<RequestPortfolio<Material>::Ptr> ports;
  RequestPortfolio<Material>::Ptr port(new RequestPortfolio<Material>());
  Material::Ptr mat = Request_();

  if (mat->quantity() > cyclus::eps_rsrc()) {
    //TODO use multiple feed commodities?
    port->AddRequest(mat, this, feed_commod);
    ports.insert(port);
  }
  return ports;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Material::Ptr MIsoEnrich::Request_() {
  double qty = std::max(0.0, feed_inv[feed_idx].capacity()
                             - feed_inv[feed_idx].quantity());
  cyclus::Composition::Ptr comp = feed_inv_comp[feed_idx];
  return cyclus::Material::CreateUntracked(qty, feed_inv_comp[feed_idx]);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr> 
    MIsoEnrich::GetMatlBids(
      cyclus::CommodMap<cyclus::Material>::type& out_requests) {
  using cyclus::Bid;
  using cyclus::BidPortfolio;
  using cyclus::CapacityConstraint;
  using cyclus::Material;
  using cyclus::Request;
  using cyclus::toolkit::MatVec;
  using cyclus::toolkit::RecordTimeSeries;

  std::set<BidPortfolio<Material>::Ptr> ports;
  
  RecordTimeSeries<double>("supply" + tails_commod, this, 
                           tails_inv.quantity());
  // TODO think about line below
  RecordTimeSeries<double>("supply" + product_commod, this,
                           feed_inv[feed_idx].quantity());

  // TODO check how reasonable this is. Apparently requests do not
  // feature the recipe, only the commodity.
  // TODO think about whether or not to implement multiple tails inventories
  if ((out_requests.count(tails_commod) > 0) 
      && (tails_inv.quantity() > 0)) {
    BidPortfolio<Material>::Ptr tails_port(new BidPortfolio<Material>());
    
    std::vector<Request<Material>*>& tails_requests = 
      out_requests[tails_commod];
    std::vector<Request<Material>*>::iterator it;
    for (it = tails_requests.begin(); it!= tails_requests.end(); it++) {
      MatVec materials = tails_inv.PopN(tails_inv.count());
      tails_inv.Push(materials);
      for (int k = 0; k < materials.size(); k++) {
        Material::Ptr m = materials[k];
        Request<Material>* req = *it;
        tails_port->AddBid(req, m, this);
      }
    }
    // TODO remove tails_inv capacity constraint and replace it
    // with multiple inventories of different compositions? 
    CapacityConstraint<Material> tails_constraint(tails_inv.quantity());
    tails_port->AddConstraint(tails_constraint);
    LOG(cyclus::LEV_INFO5, "MIsoEn") << prototype()
                                     << " adding tails capacity constraint"
                                     << " of " << tails_inv.quantity();
    ports.insert(tails_port);
  }
  
  // TODO here, one idea would be to change the if-clause and to check if 
  // any of the feed_inv's is not zero. If the current is not, then use 
  // that one, else change
  if ((out_requests.count(product_commod) > 0) 
      && (feed_inv[feed_idx].quantity() > 0)) {
    BidPortfolio<Material>::Ptr commod_port(new BidPortfolio<Material>());

    std::vector<Request<Material>*>& commod_requests = 
        out_requests[product_commod];
    std::vector<Request<Material>*>::iterator it;
    for (it = commod_requests.begin(); it != commod_requests.end(); it++) {
      Request<Material>* req = *it;
      Material::Ptr req_mat = req->target(); 
      if (ValidReq_(req_mat)) {
        Material::Ptr offer = Offer_(req_mat);
        commod_port->AddBid(req, offer, this);
      }
    }
  
    cyclus::Composition::Ptr feed_comp = feed_inv_comp[feed_idx];
    cyclus::Converter<Material>::Ptr swu_converter(
        new SwuConverter(feed_comp, tails_assay, gamma_235, 
                         use_downblending));
    cyclus::Converter<Material>::Ptr feed_converter(
        new FeedConverter(feed_comp, tails_assay, gamma_235, 
                          use_downblending));
    CapacityConstraint<Material> swu_constraint(swu_capacity, 
                                                swu_converter);
    CapacityConstraint<Material> feed_constraint(
        feed_inv[feed_idx].quantity(), feed_converter);
    commod_port->AddConstraint(swu_constraint);
    commod_port->AddConstraint(feed_constraint);

    LOG(cyclus::LEV_INFO5, "MIsoEn") << prototype() 
                                     << " adding a SWU constraint of "
                                     << swu_constraint.capacity();
    LOG(cyclus::LEV_INFO5, "MIsoEn") << prototype() 
                                     << " adding a feed constraint of "
                                     << feed_constraint.capacity();
    ports.insert(commod_port);
  }
  return ports;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Material::Ptr MIsoEnrich::Offer_(
    cyclus::Material::Ptr mat) {
  cyclus::CompMap product_comp, dummy_comp;
  double dummy_double;
  double feed_qty = feed_inv[feed_idx].quantity();
  double product_qty = mat->quantity();
  int dummy_int;
  
  enrichment_calc.SetInput(feed_inv_comp[feed_idx], MIsoAtomAssay(mat), 
                           tails_assay, feed_qty, product_qty, 
                           swu_capacity, gamma_235, use_downblending);
  enrichment_calc.EnrichmentOutput(product_comp, dummy_comp, dummy_double,
                                   dummy_double, product_qty, dummy_double,
                                   dummy_int, dummy_int);
  return cyclus::Material::CreateUntracked(
      product_qty, cyclus::Composition::CreateFromAtom(product_comp)); 
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bool MIsoEnrich::ValidReq_(const cyclus::Material::Ptr& req_mat) {
  double u_235 = MIsoAtomAssay(req_mat);
  double u_238 = MIsoAtomFrac(req_mat, IsotopeToNucID(238));

// bool u_238_present = u_238 > 0 && !cyclus::AlmostEq(u_238, 0);
  bool u_238_present = u_238 > 0;
  bool not_depleted = u_235 > tails_assay;
  bool possible_enrichment = u_235 < max_enrich;
  
  return u_238_present && not_depleted && possible_enrichment;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bool SortBids(cyclus::Bid<cyclus::Material>* i,
              cyclus::Bid<cyclus::Material>* j) {
  cyclus::Material::Ptr mat_i = i->offer();
  cyclus::Material::Ptr mat_j = j->offer();

  // TODO cycamore uses mass(U235) compared to total mass. This would also
  // include possible non-U elements that are sent directly to tails. 
  // Because of this, they should not be considered here IMO.
  return MIsoAtomAssay(mat_i) <= MIsoAtomAssay(mat_j);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrich::AdjustMatlPrefs(
    cyclus::PrefMap<cyclus::Material>::type& prefs) {
  using cyclus::Bid;
  using cyclus::Material;

  if (order_prefs == false) {
    return;
  }

  cyclus::PrefMap<cyclus::Material>::type::iterator reqit;
  // loop over all requests
  for (reqit = prefs.begin(); reqit != prefs.end(); reqit++) {
    std::vector<Bid<Material>*> bids_vector;
    std::map<Bid<Material>*, double>::iterator mit;
    // loop over all bids per request
    for (mit = reqit->second.begin(); mit != reqit->second.end(); mit++) {
      Bid<Material>* bid = mit->first;
      bids_vector.push_back(bid);
    }  // each bid
    std::sort(bids_vector.begin(), bids_vector.end(), SortBids);
    
    // The bids vector has already been sorted starting with lowest (or 
    // zero) U235 content. The following loop sets the preferences for 
    // every request with 0 U235 content to -1 such that they are ignored.
    bool u_235_mass = false;
    for (int bid_i = 0; bid_i < bids_vector.size(); bid_i++) {
      int new_pref = bid_i + 1;

      if (!u_235_mass) {
        cyclus::Material::Ptr mat = bids_vector[bid_i]->offer();
        if (MIsoAtomAssay(mat) == 0.) {
          new_pref = -1;
        } else {
          u_235_mass = true;
        }
      }
      (reqit->second)[bids_vector[bid_i]] = new_pref;
    }  // each bid
  }  // each material request

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrich::GetMatlTrades(
    const std::vector<cyclus::Trade<cyclus::Material> >& trades,
    std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                          cyclus::Material::Ptr> >& responses) {
  using cyclus::Material;
  using cyclus::Trade;

  intra_timestep_swu = 0;
  intra_timestep_feed = 0;

  std::vector<Trade<Material> >::const_iterator it;
  for (it = trades.begin(); it != trades.end(); it++) {
    double qty = it->amt;
    std::string commod_type = it->bid->request()->commodity();
    Material::Ptr response;
    
    if (commod_type == tails_commod) {
      LOG(cyclus::LEV_INFO5, "MIsoEn") << prototype() 
                                       << " just received an order for " 
                                       << it->amt << " of " 
                                       << tails_commod;
      double pop_qty = std::min(qty, tails_inv.quantity());
      response = tails_inv.Pop(pop_qty, cyclus::eps_rsrc());
    } else {
      LOG(cyclus::LEV_INFO5, "MIsoEn") << prototype() 
                                       << " just received an order for "
                                       << it->amt << " of " 
                                       << product_commod;
      response = Enrich_(it->bid->offer(), qty);
    }
    responses.push_back(std::make_pair(*it, response));
  }

  if (cyclus::IsNegative(tails_inv.quantity())) {
    std::stringstream ss;
    ss << "is being asked to provide more than its current inventory.";
    throw cyclus::ValueError(Agent::InformErrorMsg(ss.str()));
  }
  if (cyclus::IsNegative(current_swu_capacity)) {
    std::stringstream ss;
    ss << " is being asked to provide more than its SWU capacity.";
    throw cyclus::ValueError(Agent::InformErrorMsg(ss.str()));
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrich::AcceptMatlTrades(
    const std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                                cyclus::Material::Ptr> >& responses) {
  using cyclus::Material;
  using cyclus::Trade;
  // see http://stackoverflow.com/questions/5181183/boostshared-ptr-and-inheritance
  std::vector<std::pair<Trade<Material>, 
                        Material::Ptr> >::const_iterator it;
  for (it = responses.begin(); it != responses.end(); it++) {
    AddMat_(it->second);
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrich::AddMat_(cyclus::Material::Ptr mat) {
  cyclus::CompMap cm = mat->comp()->atom();
  bool non_u_elem = false;  
  
  cyclus::CompMap::const_iterator it;
  for (it = cm.begin(); it != cm.end(); it++) {
    if ((pyne::nucname::znum(it->first) != 92) && (it->second > 0)) {
      non_u_elem = true;
    }
  }
  if (non_u_elem) {
    cyclus::Warn<cyclus::VALUE_WARNING>("Non-uranium elements are sent "
                                        "directly to tails.");
  }
  AddFeedMat_(mat);

  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Material::Ptr MIsoEnrich::Enrich_(
    cyclus::Material::Ptr mat, double qty) {

  cyclus::CompMap product_comp, tails_comp;
  double feed_required, swu_required, product_qty, tails_qty;
  int n_enriching, n_stripping;

  double feed_assay = MIsoAtomAssay(feed_inv_comp[feed_idx]);
  double product_assay = MIsoAtomAssay(mat);
  
  // In the following line, the enrichment is calculated but it is not yet
  // performed!
  enrichment_calc.SetInput(feed_inv_comp[feed_idx], product_assay,
                           tails_assay, feed_inv[feed_idx].quantity(), 
                           qty, current_swu_capacity, gamma_235,
                           use_downblending);
  enrichment_calc.EnrichmentOutput(product_comp, tails_comp, feed_required,
                                   swu_required, product_qty, tails_qty,
                                   n_enriching, n_stripping);
  // Now, perform the enrichment by popping the feed and converting it to 
  // product and tails.
  cyclus::Material::Ptr pop_mat;
  try {
    if (cyclus::AlmostEq(feed_required, feed_inv[feed_idx].quantity())) {
      pop_mat = cyclus::toolkit::Squash(
          feed_inv[feed_idx].PopN(feed_inv[feed_idx].count()));
    } else {
      pop_mat = feed_inv[feed_idx].Pop(feed_required, cyclus::eps_rsrc());
    }
  } catch (cyclus::Error& e) {
    std::stringstream ss;
    ss << " tried to remove " << feed_required << " from its feed "
       << " inventory nr " << feed_idx << " holding " 
       << feed_inv[feed_idx].quantity();
    throw cyclus::ValueError(cyclus::Agent::InformErrorMsg(ss.str()));
  }
  cyclus::Material::Ptr response = pop_mat->ExtractComp(
      product_qty, cyclus::Composition::CreateFromAtom(product_comp));
  tails_inv.Push(pop_mat);

  current_swu_capacity -= swu_required;
  intra_timestep_swu += swu_required;
  intra_timestep_feed += feed_required;
  RecordEnrichment_(feed_required, swu_required, feed_idx);

  LOG(cyclus::LEV_INFO5, "MIsoEn") << prototype()
                                   << " has performed an enrichment: ";
  LOG(cyclus::LEV_INFO5, "MIsoEn") << "   * Feed Qty: " << feed_required;
  LOG(cyclus::LEV_INFO5, "MIsoEn") << "   * Feed Assay: "
                                   << feed_assay;
  LOG(cyclus::LEV_INFO5, "MIsoEn") << "   * Product Qty: " << product_qty;
  LOG(cyclus::LEV_INFO5, "MIsoEn") << "   * Product Assay: "
                                   << MIsoAtomAssay(product_comp);
  LOG(cyclus::LEV_INFO5, "MIsoEn") << "   * Tails Qty: " << tails_qty;
  LOG(cyclus::LEV_INFO5, "MIsoEn") << "   * Tails Assay: "
                                   << MIsoAtomAssay(tails_comp);
  LOG(cyclus::LEV_INFO5, "MIsoEn") << "   * SWU: " << swu_required;
  LOG(cyclus::LEV_INFO5, "MIsoEn") << "   * Current SWU capacity: " 
                                   << current_swu_capacity;

  return response;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrich::RecordEnrichment_(double feed_qty, double swu,
                                           int feed_inv_idx) {
  LOG(cyclus::LEV_DEBUG1, "MIsoEn") << prototype()
                                    << " has enriched a material:";
  LOG(cyclus::LEV_DEBUG1, "MIsoEn") << "  * Amount: " << feed_qty;
  LOG(cyclus::LEV_DEBUG1, "MIsoEn") << "  *    SWU: " << swu;

  cyclus::Context* ctx = cyclus::Agent::context();
  ctx->NewDatum("MIsoEnrichments")
     ->AddVal("AgentId", id())
     ->AddVal("Time", ctx->time())
     ->AddVal("feed_qty", feed_qty)
     ->AddVal("feed_inv_idx", feed_inv_idx)
     ->AddVal("SWU", swu)
     ->Record();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void MIsoEnrich::RecordPosition() {
  std::string specification = this->spec();
  context()->NewDatum("AgentPosition")
           ->AddVal("Spec", specification)
           ->AddVal("Prototype", this->prototype())
           ->AddVal("AgentId", id())
           ->AddVal("Latitude", latitude)
           ->AddVal("Longitude", longitude)
           ->Record();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// WARNING! Do not change the following this function!!! This enables your
// archetype to be dynamically loaded and any alterations will cause your
// archetype to fail.
extern "C" cyclus::Agent* ConstructMIsoEnrich(cyclus::Context* ctx) {
  return new MIsoEnrich(ctx);
}

}  // namespace misoenrichment
