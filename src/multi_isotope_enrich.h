#ifndef MULTIISOTOPEENRICHMENT_SRC_MULTI_ISOTOPE_ENRICH_H_
#define MULTIISOTOPEENRICHMENT_SRC_MULTI_ISOTOPE_ENRICH_H_

#include <string>
#include <vector>

#include "cyclus.h"

namespace multiisotopeenrichment {

/// @class MultiIsotopeEnrich
///
/// This Facility is intended
/// as a skeleton to guide the implementation of new Facility
/// agents.
/// The MultiIsotopeEnrich class inherits from the Facility class and is
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
class MultiIsotopeEnrich : public cyclus::Facility,
                           public cyclus::toolkit::Position {
 public:
  /// Constructor for MultiIsotopeEnrich Class
  /// @param ctx the cyclus context for access to simulation-wide parameters
  explicit MultiIsotopeEnrich(cyclus::Context* ctx);

  /// The Prime Directive
  /// @warning The Prime Directive must have a space before it! (A fix will be
  /// in 2.0 ^TM)

  #pragma cyclus

  #pragma cyclus note {"doc": "A stub facility is provided as a skeleton " \
                              "for the design of new facility agents."}

  virtual std::string str();
  virtual void Build(cyclus::Agent* parent);
  
  virtual void Tick();
  virtual void Tock();
  
  virtual std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr>
      GetMatlRequests();

  virtual void AdjustMatlPrefs(cyclus::PrefMap<cyclus::Material>::type& prefs);

  virtual void AcceptMatlTrades(
      const std::vector< std::pair<cyclus::Trade<cyclus::Material>,
      cyclus::Material::Ptr> >& responses);

  virtual std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr>
    GetMatlBids(cyclus::CommodMap<cyclus::Material>::type&
    commod_requests);

  virtual void GetMatlTrades(
    const std::vector< cyclus::Trade<cyclus::Material> >& trades,
    std::vector<std::pair<cyclus::Trade<cyclus::Material>,
    cyclus::Material::Ptr> >& responses);

  bool ValidReq(const cyclus::Material::Ptr mat);

 private:
  void AddMat_(cyclus::Material::Ptr mat);

  cyclus::Material::Ptr Request_();
  
  cyclus::Material::Ptr Offer_(cyclus::Material::Ptr req);

  cyclus::Material::Ptr Enrich_(cyclus::Material::Ptr mat, double qty);

  double FeedAssay();

  ///  @brief records and enrichment with the cyclus::Recorder
  void RecordEnrichment_(double natural_u, double swu);

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
  double max_swu_capacity;
  double current_swu_capacity;

  double intra_timestep_swu;
  double intra_timestep_feed;
  
  std::vector<cyclus::toolkit::ResBuf<cyclus::Material>> feed_inv;
  std::vector<cyclus::toolkit::ResBuf<cyclus::Material>> tails_inv;

  std::vector<cyclus::Composition::Ptr> feed_inv_comp;
  int current_feed_inv;
  
  
};

}  // namespace multiisotopeenrichment

#endif  // CYCLUS_MULTIISOTOPEENRICHMENT_MULTI_ISOTOPE_ENRICH_H_
