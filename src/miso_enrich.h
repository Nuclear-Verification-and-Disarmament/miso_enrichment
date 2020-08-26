#ifndef MISOENRICHMENT_SRC_MISO_ENRICH_H_
#define MISOENRICHMENT_SRC_MISO_ENRICH_H_

#include <string>
#include <vector>

#include "cyclus.h"

#include "enrichment_calculator.h"
#include "flexible_input.cc"
#include "miso_helper.h"

namespace misoenrichment {

// TODO think about passing a pointer to MIsoEnrich's EnrichmentCalculator
// object? Think about how this works with possibly uninitialised gamma_235
class SwuConverter : public cyclus::Converter<cyclus::Material> {
 public:
  SwuConverter(cyclus::Composition::Ptr feed_comp, double tails_assay,
               double gamma_235) 
      : feed_comp_(feed_comp), gamma_235_(gamma_235),
        tails_assay_(tails_assay) {}
  
  virtual ~SwuConverter() {}

  virtual double convert(
      cyclus::Material::Ptr m, cyclus::Arc const * a = NULL,
      cyclus::ExchangeTranslationContext<cyclus::Material> 
          const * ctx = NULL) const {
    EnrichmentCalculator e;
    
    double product_qty = m->quantity();
    double product_assay = MIsoAtomAssay(m);
    e.SetInput(feed_comp_, product_assay, tails_assay_, 1e299, product_qty,
               1e299, gamma_235_);
    double swu_used = e.SwuUsed();
    
    return swu_used;
  }

  virtual bool operator==(Converter& other) const {
    SwuConverter* cast = dynamic_cast<SwuConverter*>(&other);
    
    bool cast_not_null = cast != NULL;
    bool feed_eq = cyclus::compmath::AlmostEq(feed_comp_->atom(), 
                                              cast->feed_comp_->atom(),
                                              kEpsCompMap);
    bool tails_eq = tails_assay_ == cast->tails_assay_;

    return cast != NULL && feed_eq && tails_eq;
  }

 private:
  cyclus::Composition::Ptr feed_comp_;
  double gamma_235_;
  double tails_assay_;
};

// TODO think about passing a pointer to MIsoEnrich's EnrichmentCalculator
// object? Think about how this works with possibly uninitialised gamma_235
class FeedConverter : public cyclus::Converter<cyclus::Material> {
 public:
  FeedConverter(cyclus::Composition::Ptr feed_comp, double tails_assay,
                double gamma_235)
      : feed_comp_(feed_comp), gamma_235_(gamma_235), 
        tails_assay_(tails_assay) {}

  virtual ~FeedConverter() {}

  virtual double convert(
      cyclus::Material::Ptr m, cyclus::Arc const * a = NULL,
      cyclus::ExchangeTranslationContext<cyclus::Material> 
          const * ctx = NULL) const {
    
    double product_qty = m->quantity();
    double product_assay = MIsoAtomAssay(m);
    EnrichmentCalculator e;
    e.SetInput(feed_comp_, product_assay, tails_assay_, 1e299, product_qty,
               1e299, gamma_235_);
    double feed_used = e.FeedUsed();
    
    cyclus::toolkit::MatQuery mq(m);
    std::vector<int> isotopes;
    IsotopesNucID(isotopes);
    std::set<int> nucs(isotopes.begin(), isotopes.end());
    double feed_uranium_frac = mq.atom_frac(nucs);

    return feed_used / feed_uranium_frac;
  }
  
  virtual bool operator==(Converter& other) const {
    FeedConverter* cast = dynamic_cast<FeedConverter*>(&other);
    
    bool cast_not_null = cast != NULL;
    bool feed_eq = cyclus::compmath::AlmostEq(feed_comp_->atom(), 
                                              cast->feed_comp_->atom(),
                                              kEpsCompMap);
    bool tails_eq = tails_assay_ == cast->tails_assay_;

    return cast != NULL && feed_eq && tails_eq;
  }

 private:
  cyclus::Composition::Ptr feed_comp_;
  double gamma_235_;
  double tails_assay_;
};

/// @class MIsoEnrich
///
/// @section intro 
///
/// @section agentparams 
///
/// @section optionalparams
///
/// @section detailed 
class MIsoEnrich : public cyclus::Facility,
                   public cyclus::toolkit::Position {
 public:
  /// Constructor for MIsoEnrich Class
  /// @param ctx the cyclus context for access to simulation-wide parameters
  explicit MIsoEnrich(cyclus::Context* ctx);

  ~MIsoEnrich();
  
  friend class MIsoEnrichTest;

  #pragma cyclus
  #pragma cyclus note {"doc": "A stub facility is provided as a skeleton " \
                              "for the design of new facility agents."}
  
  inline void EnterNotify() {
    cyclus::Facility::EnterNotify();
  };
  void Build(cyclus::Agent* parent);
  void Tick();
  void Tock();
  void AdjustMatlPrefs(cyclus::PrefMap<cyclus::Material>::type& prefs);
  void AcceptMatlTrades(
      const std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                                  cyclus::Material::Ptr> >& responses);
  void GetMatlTrades(
      const std::vector<cyclus::Trade<cyclus::Material> >& trades,
      std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                            cyclus::Material::Ptr> >& responses);


  std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr>
      GetMatlBids(cyclus::CommodMap<cyclus::Material>::type& 
                  commod_requests);
  std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr> 
      GetMatlRequests();
  std::string str();

 private:
  void AddMat_(cyclus::Material::Ptr mat);
  
  void AddFeedMat_(cyclus::Material::Ptr mat);

  cyclus::Material::Ptr Request_();

  // The Offer function only considers U235 content that needs to be 
  // achieved and it ignores the minor isotopes. This has the advantage 
  // that the evolution of minor isotopes does not need to be taken into 
  // account when performing requests to a MIsoEnrich facility.
  cyclus::Material::Ptr Offer_(cyclus::Material::Ptr req);

  cyclus::Material::Ptr Enrich_(cyclus::Material::Ptr mat, double qty);

  bool ValidReq_(const cyclus::Material::Ptr& mat);
  
  ///  @brief records and enrichment with the cyclus::Recorder
  void RecordEnrichment_(double feed_qty, double swu, int feed_inv_idx);

  /// Records an agent's latitude and longitude to the output db
  void RecordPosition();

  #pragma cyclus var { \
    "tooltip": "feed commodity", \
    "doc": "feed commodity that the enrichment facility accepts",	\
    "uilabel": "Feed Commodity", \
    "uitype": "incommodity" \
  }
  std::string feed_commod;

  #pragma cyclus var { \
    "tooltip": "feed recipe",	\
    "doc": "recipe for enrichment facility feed commodity", \
    "uilabel": "Feed Recipe", \
    "uitype": "inrecipe" \
  }
  std::string feed_recipe;

  #pragma cyclus var { \
    "tooltip": "product commodity",	\
    "doc": "product commodity that the enrichment facility generates", \
    "uilabel": "Product Commodity", \
    "uitype": "outcommodity" \
  }
  std::string product_commod;

  #pragma cyclus var { \
    "tooltip": "tails commodity",	\
    "doc": "tails commodity supplied by enrichment facility", \
    "uilabel": "Tails Commodity", \
    "uitype": "outcommodity" \
  }
  std::string tails_commod;

  #pragma cyclus var { \
    "default": 0.003, "tooltip": "tails assay",	\
    "uilabel": "Tails Assay", \
    "uitype": "range", \
    "doc": "tails assay from the enrichment process", \
  }
  double tails_assay;

  #pragma cyclus var { \
    "default": 0, "tooltip": "initial uranium reserves (kg)",	\
    "uilabel": "Initial Feed Inventory", \
    "doc": "amount of natural uranium stored at the enrichment " \
    "facility at the beginning of the simulation (kg)" \
  }
  double initial_feed;

  #pragma cyclus var { \
    "default": 1e299, "tooltip": "max inventory of feed material (kg)", \
    "uilabel": "Maximum Feed Inventory", \
    "uitype": "range", \
    "range": [0.0, 1e299], \
    "doc": "maximum total inventory of natural uranium in "	\
           "the enrichment facility (kg)" \
  }
  double max_feed_inventory;

  #pragma cyclus var { \
    "default": 1.0,	\
    "tooltip": "maximum allowed enrichment fraction", \
    "doc": "maximum allowed weight fraction of U235 in product", \
    "uilabel": "Maximum Allowed Enrichment", \
    "uitype": "range", \
    "range": [0.0,1.0], \
    "schema": '<optional>' \
        '          <element name="max_enrich">' \
        '              <data type="double">' \
        '                  <param name="minInclusive">0</param>' \
        '                  <param name="maxInclusive">1</param>' \
        '              </data>'	\
        '          </element>' \
        '      </optional>'	\
  }
  double max_enrich;

  #pragma cyclus var { \
    "default": 1,	\
    "userlevel": 10, \
    "tooltip": "Rank Material Requests by U235 Content", \
    "uilabel": "Prefer feed with higher U235 content", \
    "doc": "turn on preference ordering for input material " \
           "so that EF chooses higher U235 content first" \
  }
  bool order_prefs;

  #pragma cyclus var { \
    "default": 1.4, \
    "tooltip": "Separation factor U235", \
    "uilabel": "Separation factor for U235", \
    "doc": "overall stage separation factor for U235" \
  }
  double gamma_235;

  #pragma cyclus var { \
    "default": 1e299,	\
    "tooltip": "SWU capacity (kgSWU/month)", \
    "uilabel": "SWU Capacity", \
    "uitype": "range", \
    "range": [0.0, 1e299], \
    "doc": "separative work unit (SWU) capacity of enrichment "	\
           "facility (kgSWU/timestep) " \
  }
  double swu_capacity;
  double current_swu_capacity;

  double intra_timestep_swu;
  double intra_timestep_feed;

  EnrichmentCalculator enrichment_calc;
  
  // TODO think about how to include these variables in preprocessor
  //#pragma cyclus var {}
  std::vector<cyclus::toolkit::ResBuf<cyclus::Material> > feed_inv;
  //#pragma cyclus var {}
  std::vector<cyclus::Composition::Ptr> feed_inv_comp;

  int feed_idx;
  
  #pragma cyclus var {}
  cyclus::toolkit::ResBuf<cyclus::Material> tails_inv;

  #pragma cyclus var { \
    "default": 0.0, \
    "uilabel": "Geographical latitude in degrees as a double", \
    "doc": "Latitude of the agent's geographical position. The value " \
           "should be expressed in degrees as a double." \
  } 
  double latitude;
  
  #pragma cyclus var { \
    "default": 0.0, \
    "uilabel": "Geographical longitude in degrees as a double", \
    "doc": "Longitude of the agent's geographical position. The value " \
           "should be expressed in degrees as a double." \
  } 
  double longitude; 
  
  cyclus::toolkit::Position coordinates;
  
  #pragma cyclus var {  \
    "default": [],  \
    "tooltip": "SWU capacity change times in timesteps from beginning " \
               "of deployment",  \
    "uilabel": "SWU capacity change times",  \
    "uitype": "oneormore",  \
    "doc": "list of timesteps where the SWU is changed. The first " \
           "timestep has to be 0 as it sets the initial value and all " \
           "timesteps are measured from the moment of deployment of the " \
           "facility, not from the start of the simulation."  \
  }
  std::vector<int> swu_capacity_times;

  #pragma cyclus var {  \
    "default": [],  \
    "tooltip": "SWU capacity list (kg SWU / month)",  \
    "uilabel": "SWU Capacity List",  \
    "uitype": "oneormore",  \
    "doc": "list of separative work unit (SWU) capacity of enrichment " \
           "facility (kg SWU / month)"  \
  }
  std::vector<double> swu_capacity_vals;
  FlexibleInput<double> swu_flexible;

};

}  // namespace misoenrichment

#endif  // MISOENRICHMENT_MISO_ENRICH_H_
