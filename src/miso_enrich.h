#ifndef MISOENRICHMENT_SRC_MISO_ENRICH_H_
#define MISOENRICHMENT_SRC_MISO_ENRICH_H_

#include <string>
#include <vector>

#include "cyclus.h"

#include "enrichment_calculator.h"
#include "miso_helper.h"

namespace misoenrichment {

class SwuConverter : public cyclus::Converter<cyclus::Material> {
 public:
  SwuConverter(cyclus::Composition::Ptr feed_comp, double tails_assay) 
      : feed_comp_(feed_comp),
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
               1e299);
    double swu_used = e.SwuUsed();
    
    return swu_used;
  }

  virtual bool operator==(Converter& other) const {
    SwuConverter* cast = dynamic_cast<SwuConverter*>(&other);
    
    // TODO make a class for the whole module with limits?
    const double kEpsComp = 1e-10;

    bool cast_not_null = cast != NULL;
    bool feed_eq = cyclus::compmath::AlmostEq(feed_comp_->atom(), 
                                              cast->feed_comp_->atom(),
                                              kEpsComp);
    bool tails_eq = tails_assay_ == cast->tails_assay_;

    return cast != NULL && feed_eq && tails_eq;
  }

 private:
  cyclus::Composition::Ptr feed_comp_;
  double tails_assay_;
};

class FeedConverter : public cyclus::Converter<cyclus::Material> {
 public:
  FeedConverter(cyclus::Composition::Ptr feed_comp, double tails_assay)
      : feed_comp_(feed_comp),
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
               1e299);
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
    
    // TODO make a class for the whole module with limits?
    const double kEpsComp = 1e-10;

    bool cast_not_null = cast != NULL;
    bool feed_eq = cyclus::compmath::AlmostEq(feed_comp_->atom(), 
                                              cast->feed_comp_->atom(),
                                              kEpsComp);
    bool tails_eq = tails_assay_ == cast->tails_assay_;

    return cast != NULL && feed_eq && tails_eq;
  }

 private:
  cyclus::Composition::Ptr feed_comp_;
  double tails_assay_;
};
/// @class MIsoEnrich
///
/// This Facility is intended
/// as a skeleton to guide the implementation of new Facility
/// agents.
/// The MIsoEnrich class inherits from the Facility class and is
/// dynamically loaded by the Agent class when requested.
///
/// @section intro Introduction
/// Place an introduction to the agent here.
///
/// @section agentparams Agent Parameters
/// Place a description of the required input parameters which define the
/// agent implementation.
///
/// @section optionalparams Optional Parameters
/// Place a description of the optional input parameters to define the
/// agent implementation.
///
/// @section detailed Detailed Behavior
/// Place a description of the detailed behavior of the agent. Consider
/// describing the behavior at the tick and tock as well as the behavior
/// upon sending and receiving materials and messages.
class MIsoEnrich : public cyclus::Facility,
                           public cyclus::toolkit::Position {
 public:
  /// Constructor for MIsoEnrich Class
  /// @param ctx the cyclus context for access to simulation-wide parameters
  explicit MIsoEnrich(cyclus::Context* ctx);

  /// The Prime Directive
  /// @warning The Prime Directive must have a space before it! (A fix will be
  /// in 2.0 ^TM)
  #pragma cyclus

  #pragma cyclus note {"doc": "A stub facility is provided as a skeleton " \
                              "for the design of new facility agents."}
  
  virtual ~MIsoEnrich();

  virtual std::string str();
  virtual void Build(cyclus::Agent* parent);
  
  virtual void Tick();
  virtual void Tock();
  
  virtual std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr>
      GetMatlRequests();

  virtual void AdjustMatlPrefs(
      cyclus::PrefMap<cyclus::Material>::type& prefs);

  virtual void AcceptMatlTrades(
      const std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                                  cyclus::Material::Ptr> >& responses);

  virtual std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr>
      GetMatlBids(cyclus::CommodMap<cyclus::Material>::type& 
                  commod_requests);

  virtual void GetMatlTrades(
    const std::vector<cyclus::Trade<cyclus::Material> >& trades,
    std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                          cyclus::Material::Ptr> >& responses);

  bool ValidReq(const cyclus::Material::Ptr mat);

 private:
  void AddMat_(cyclus::Material::Ptr mat);

  cyclus::Material::Ptr Request_();

  // The Offer function only considers U235 content that needs to be 
  // achieved and it ignores the minor isotopes. This has the advantage 
  // that the evolution of minor isotopes does not need to be taken into 
  // account when performing requests to a MIsoEnrich facility.
  cyclus::Material::Ptr Offer_(cyclus::Material::Ptr req);

  cyclus::Material::Ptr Enrich_(cyclus::Material::Ptr mat, double qty);

  bool ValidReq_(const cyclus::Material::Ptr mat);
  
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
    "uilabel": "Feed Recipe",                                   \
    "uitype": "inrecipe" \
  }
  std::string feed_recipe;

  #pragma cyclus var { \
    "tooltip": "product commodity",					\
    "doc": "product commodity that the enrichment facility generates",	 \
    "uilabel": "Product Commodity",                                     \
    "uitype": "outcommodity" \
  }
  std::string product_commod;

  #pragma cyclus var {							\
    "tooltip": "tails commodity",					\
    "doc": "tails commodity supplied by enrichment facility",		\
    "uilabel": "Tails Commodity",                                   \
    "uitype": "outcommodity" \
  }
  std::string tails_commod;

  #pragma cyclus var {							\
    "default": 0.003, "tooltip": "tails assay",				\
    "uilabel": "Tails Assay",                             \
    "uitype": "range",                               \
    "doc": "tails assay from the enrichment process",       \
  }
  double tails_assay;

  #pragma cyclus var {							\
    "default": 0, "tooltip": "initial uranium reserves (kg)",		\
    "uilabel": "Initial Feed Inventory",				\
    "doc": "amount of natural uranium stored at the enrichment "	\
    "facility at the beginning of the simulation (kg)"			\
  }
  double initial_feed;

  #pragma cyclus var {							\
    "default": 1e299, "tooltip": "max inventory of feed material (kg)", \
    "uilabel": "Maximum Feed Inventory", \
    "uitype": "range", \
    "range": [0.0, 1e299], \
    "doc": "maximum total inventory of natural uranium in "		\
           "the enrichment facility (kg)"     \
  }
  double max_feed_inventory;

  #pragma cyclus var { \
    "default": 1.0,						\
    "tooltip": "maximum allowed enrichment fraction",		\
    "doc": "maximum allowed weight fraction of U235 in product", \
    "uilabel": "Maximum Allowed Enrichment", \
    "uitype": "range", \
    "range": [0.0,1.0], \
    "schema": '<optional>'				     	   \
        '          <element name="max_enrich">'			   \
        '              <data type="double">'			   \
        '                  <param name="minInclusive">0</param>'   \
        '                  <param name="maxInclusive">1</param>'   \
        '              </data>'					   \
        '          </element>'					   \
        '      </optional>'					   \
  }
  double max_enrich;

  #pragma cyclus var { \
    "default": 1,		       \
    "userlevel": 10,							\
    "tooltip": "Rank Material Requests by U235 Content",		\
    "uilabel": "Prefer feed with higher U235 content", \
    "doc": "turn on preference ordering for input material "		\
           "so that EF chooses higher U235 content first" \
  }
  bool order_prefs;

  #pragma cyclus var {						       \
    "default": 1e299,						       \
    "tooltip": "SWU capacity (kgSWU/month)",			       \
    "uilabel": "SWU Capacity",                                         \
    "uitype": "range",                                                  \
    "range": [0.0, 1e299],                                               \
    "doc": "separative work unit (SWU) capacity of enrichment "		\
           "facility (kgSWU/timestep) "                                     \
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
  //#pragma cyclus var {}
  int feed_idx;
  
  //#pragma cyclus var {}
  cyclus::toolkit::ResBuf<cyclus::Material> tails_inv;

  #pragma cyclus var { \
    "default": 0.0, \
    "uilabel": "Geographical latitude in degrees as a double", \
    "doc": "Latitude of the agent's geographical position. The value " \
           " should be expressed in degrees as a double." \
  } 
  double latitude;
  
  #pragma cyclus var { \
    "default": 0.0, \
    "uilabel": "Geographical longitude in degrees as a double", \
    "doc": "Longitude of the agent's geographical position. The value " \
           " should be expressed in degrees as a double." \
  } 
  double longitude; 
  
  cyclus::toolkit::Position coordinates;
};

}  // namespace misoenrichment

#endif  // MISOENRICHMENT_MISO_ENRICH_H_
