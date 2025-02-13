#include "var_recipe_source.h"

#include <algorithm>  // std::min
#include <cstring>  // std::memcpy
#include <sstream>
#include <limits>

#include "pyne.h"  // pyne::nucname::isnuclide

namespace misoenrichment {

VarRecipeSource::VarRecipeSource(cyclus::Context* ctx)
    : cyclus::Facility(ctx),
      out_commod(""),
      var_out_recipe_(
          std::pair<std::string,
                    std::map<int,
                             std::pair<std::string,
                                       std::vector<double> > > >()),
      normalisation_nuc_id_(0),
      inventory_size(1e299),
      throughput_times(std::vector<int>({-1})),
      throughput_vals(std::vector<double>({1e299})),
      flexible_throughput(FlexibleInput<double>()),
      current_throughput(0.),
      latitude(0.0),
      longitude(0.0),
      coordinates(latitude, longitude) {}

VarRecipeSource::~VarRecipeSource() {}

void VarRecipeSource::InitFrom(VarRecipeSource* m) {
  #pragma cyclus impl initfromcopy misoenrichment::VarRecipeSource
  cyclus::toolkit::CommodityProducer::Copy(m);
  RecordPosition_();
}

void VarRecipeSource::InitFrom(cyclus::QueryableBackend* b) {
  #pragma cyclus impl initfromdb misoenrichment::VarRecipeSource
  cyclus::toolkit::CommodityProducer::Add(
      cyclus::toolkit::Commodity(out_commod),
      cyclus::toolkit::CommodInfo(current_throughput, current_throughput));
  RecordPosition_();
}

std::string VarRecipeSource::str() {
  std::stringstream ss;
  std::string ans;
  if (cyclus::toolkit::CommodityProducer::Produces(
          cyclus::toolkit::Commodity(out_commod))) {
    ans = "yes";
  } else {
    ans = "no";
  }
  ss << cyclus::Facility::str()
     << " supplies commodity '" << out_commod
     << "' with a varying recipe (";

  for (auto const& [nuc_id, rng_properties] : var_out_recipe_.second) {
    ss << nuc_id << ": "
       << rng_properties.first << "(";
    for (double param : rng_properties.second) {
      ss << param << ", ";
    }
    ss.seekp(-2, ss.cur);  // Remove the last, superfluous ', '.
    ss << "); ";
  }
  ss.seekp(-2, ss.cur);
  ss << " at a current throughput of " << current_throughput
     << " kg per time step. commod producer members: "
     << " produces " << out_commod << "?: " << ans
     << " throughput: " << cyclus::toolkit::CommodityProducer::Capacity(out_commod)
     << " cost: " << cyclus::toolkit::CommodityProducer::Cost(out_commod);
  return ss.str();
}

void VarRecipeSource::EnterNotify() {
  namespace tk = cyclus::toolkit;

  cyclus::Facility::EnterNotify();

  CheckVarOutRecipe_();

  if (throughput_times[0]==-1) {
    flexible_throughput = FlexibleInput<double>(this, throughput_vals);
  } else {
    flexible_throughput = FlexibleInput<double>(this, throughput_vals,
                                                throughput_times);
  }
  current_throughput = throughput_vals[0];
  tk::CommodityProducer::SetCapacity(tk::Commodity(out_commod),
                                     current_throughput);
  tk::CommodityProducer::SetCost(tk::Commodity(out_commod), current_throughput);

  LOG(cyclus::LEV_DEBUG2, "FlxSrc") << "VarRecipeSource entering the "
                                    << "simulation: ";
  LOG(cyclus::LEV_DEBUG2, "FlxSrc") << str();
  RecordPosition_();
}

void VarRecipeSource::Tick() {
  namespace tk = cyclus::toolkit;

  // For an unknown reason, 'UpdateValue' has to be called with a copy of
  // the 'this' pointer. When directly using 'this', the address passed to
  // the function is increased by 8 bits resulting later on in a
  // segmentation fault.
  // TODO The problem described above has not been checked in VarRecipeSource
  // but it was the case in MIsoEnrichment. Maybe check if this has changed
  // (for whatever reason) here in VarRecipeSource?
  cyclus::Agent* copy_ptr;
  cyclus::Agent* source_ptr = this;
  std::memcpy((void*) &copy_ptr, (void*) &source_ptr, sizeof(cyclus::Agent*));
  current_throughput = flexible_throughput.UpdateValue(copy_ptr);
  tk::CommodityProducer::SetCapacity(tk::Commodity(out_commod),
                                     current_throughput);
  tk::CommodityProducer::SetCost(tk::Commodity(out_commod), current_throughput);

  current_composition_ = CreateRandomComposition_();
}

void VarRecipeSource::Tock() {}

std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr>
VarRecipeSource::GetMatlBids(
    cyclus::CommodMap<cyclus::Material>::type& commod_requests) {
  using cyclus::Bid;
  using cyclus::BidPortfolio;
  using cyclus::CapacityConstraint;
  using cyclus::Material;
  using cyclus::Request;

  double max_qty = std::min(current_throughput, inventory_size);
  cyclus::toolkit::RecordTimeSeries<double>("supply"+out_commod, this,
                                            max_qty);
  LOG(cyclus::LEV_INFO3, "FlxSrc") << prototype()
      << " is bidding up to " << max_qty << " kg of " << out_commod;

  std::set<BidPortfolio<Material>::Ptr> ports;
  if (max_qty < cyclus::eps()) {
    return ports;
  } else if (commod_requests.count(out_commod) == 0) {
    return ports;
  }

  BidPortfolio<Material>::Ptr port(new BidPortfolio<Material>());
  std::vector<Request<Material>*>& requests = commod_requests[out_commod];
  std::vector<Request<Material>*>::iterator it;
  for (it = requests.begin(); it != requests.end(); ++it) {
    Request<Material>* req = *it;
    Material::Ptr target = req->target();
    double qty = std::min(target->quantity(), max_qty);
    Material::Ptr m = Material::CreateUntracked(qty, current_composition_);
    port->AddBid(req, m, this);
  }

  CapacityConstraint<Material> cc(max_qty);
  port->AddConstraint(cc);
  ports.insert(port);
  return ports;
}

void VarRecipeSource::GetMatlTrades(
    const std::vector<cyclus::Trade<cyclus::Material> >& trades,
    std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                          cyclus::Material::Ptr> >& responses) {
  std::vector<cyclus::Trade<cyclus::Material> >::const_iterator it;
  for (it = trades.begin(); it != trades.end(); ++it) {
    double qty = it->amt;
    inventory_size -= qty;

    cyclus::Material::Ptr response;
    response = cyclus::Material::Create(this, qty, current_composition_);
    responses.push_back(std::make_pair(*it, response));
    LOG(cyclus::LEV_INFO5, "FlxSrc") << prototype() << " sent an order"
                                     << " for " << qty << " of " << out_commod;
  }
}

void VarRecipeSource::CheckVarOutRecipe_() {
  if (var_out_recipe_.first != "mass" && var_out_recipe_.first != "atom") {
    throw cyclus::ValueError("'mass_or_atom' must be 'mass' or 'atom'");
  }

  int n_normalisation_distributions = 0;
  for (auto const& [nuc_id, rng_properties] : var_out_recipe_.second) {
    // Check if nuc_ids are correct
    if (!pyne::nucname::isnuclide(nuc_id)) {
      std::stringstream ss;
      ss << "Nuclide id '" << nuc_id << "' is not a valid nuclide!";
      throw cyclus::ValueError(ss.str());
    }
    // Check if all distributions are valid.
    if (rng_properties.first == "normalisation") {
      n_normalisation_distributions++;
      normalisation_nuc_id_ = nuc_id;
    } else if (
        rng_properties.first != "uniform"
        && rng_properties.first != "normal") {
      throw cyclus::ValueError(
          "Distributions must be 'uniform' or 'normal'. Moreover, exactly one "
          "distribution must be 'normalisation'."
      );
    }
  }
  // Ensure that there is exactly one nuclide used for normalisation.
  if (n_normalisation_distributions != 1) {
    throw cyclus::ValueError(
        "Exactly one distribution must be 'normalisation'");
  }
}

cyclus::Composition::Ptr VarRecipeSource::CreateRandomComposition_() {
  cyclus::CompMap out_recipe_comp_map;
  double normalisation = 1.;
  // Add all nuclides to the CompMap except for the normalising nuclide.
  for (auto const& [nuc_id, rng_properties] : var_out_recipe_.second) {
    double fraction;
    if (rng_properties.first == "uniform") {
      fraction = context()->random_uniform_real(rng_properties.second[0],
                                                rng_properties.second[1]);
    } else if (rng_properties.first == "normal") {
      fraction = context()->random_normal_real(rng_properties.second[0],
                                               rng_properties.second[1],
                                               rng_properties.second[2],
                                               rng_properties.second[3]);
    } else {
      // Distribution must be 'normalisation' (the correct form of the input
      // was checked in CheckVarOutRecipe_).
      continue;
    }
    out_recipe_comp_map[nuc_id] = fraction;
    normalisation -= fraction;
  }
  out_recipe_comp_map[normalisation_nuc_id_] = normalisation;

  if (var_out_recipe_.first == "mass") {
    return cyclus::Composition::CreateFromMass(out_recipe_comp_map);
  } else {
    // No need for else if, as we have checked the input in CheckVarOutRecipe_.
    return cyclus::Composition::CreateFromAtom(out_recipe_comp_map);
  }
}

void VarRecipeSource::RecordPosition_() {
  std::string specification = this->spec();
  context()
      ->NewDatum("AgentPosition")
      ->AddVal("Spec", specification)
      ->AddVal("Prototype", this->prototype())
      ->AddVal("AgentId", id())
      ->AddVal("Latitude", latitude)
      ->AddVal("Longitude", longitude)
      ->Record();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// WARNING! Do not change the following this function!!! This enables your
// archetype to be dynamically loaded and any alterations will cause your
// archetype to fail.
extern "C" cyclus::Agent* ConstructVarRecipeSource(cyclus::Context* ctx) {
  return new VarRecipeSource(ctx);
}

}  // namespace misoenrichment
