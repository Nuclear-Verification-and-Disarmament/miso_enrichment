#ifndef MISOENRICHMENT_SRC_GPR_REACTOR_H_
#define MISOENRICHMENT_SRC_GPR_REACTOR_H_

#include <string>
#include <vector>

#include "cyclus.h"

// Future changes relating to the implementation of Antonio's GPRs are marked
// with the following comment:
// TODO ANTONIO GPR
//
// TODO list:
// - check line 49 in .cc file. Ensure that `unique_out_commods.empty()`
//   evaluates to `true`, else, the set gets not filled!
// - no implementation of `Reactor::fuel_pref(cyclus::Material::Ptr)` because
//   it is not used in the reactor class.
// - think about the weird (?) decommissioning behaviour of the cycamore
//   archetype and whether or not I should use it as well
// - check if GPRs use mass or molar fractions?

namespace misoenrichment {

class GprReactor : public cyclus::Facility, public cyclus::toolkit::Position  {
 public:
  /// Constructor for GprReactor Class
  /// @param ctx the cyclus context for access to simulation-wide parameters
  explicit GprReactor(cyclus::Context* ctx);

  ~GprReactor();

  friend class GprReactorTest;

  /// The Prime Directive
  /// Generates code that handles all input file reading and restart operations
  /// (e.g., reading from the database, instantiating a new object, etc.).
  /// @warning The Prime Directive must have a space before it!

  #pragma cyclus

  #pragma cyclus note {"doc": "A reactor facility that calculates its spent " \
                              "fuel composition using Gaussian process " \
                              "regression."}

  bool CheckDecommissionCondition();
  std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr> GetMatlBids(
      cyclus::CommodMap<cyclus::Material>::type& commod_requests);
  std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr> GetMatlRequests();
  std::string str();
  void AcceptMatlTrades(
      const std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                                  cyclus::Material::Ptr> >& responses);
  void EnterNotify();
  void GetMatlTrades(
      const std::vector<cyclus::Trade<cyclus::Material> >& trades,
      std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                            cyclus::Material::Ptr> >& responses);
  void Tick();
  void Tock();

 private:
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Fuel commodities, recipes and preferences
  #pragma cyclus var { \
    "uitype": ["oneormore", "incommodity"], \
    "doc": "An ordered list of the input commodities (fresh fuel)." \
  }
  std::vector<std::string> in_commods;

  #pragma cyclus var { \
    "uitype": ["oneormore", "outcommodity"], \
    "doc": "An ordered list of the (spent) output commodities." \
  }
  std::vector<std::string> out_commods;

  #pragma cyclus var { \
    "uitype": ["oneormore", "inrecipe"], \
    "doc": "An ordered list of the input (fresh fuel) recipes." \
  }
  std::vector<std::string> in_recipes;

  #pragma cyclus var { \
    "default": [], \
    "doc": "The preference for each input fuel type, in the same order as " \
           "the input commodities. If no preferences are specified, then the " \
           "same default value is used for all fuel requests." \
  }
  std::vector<double> fuel_prefs;

  std::set<std::string> unique_out_commods;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Core specifics
  #pragma cyclus var { \
    "default": 3, \
    "doc": "Number of assemblies in a full core. This value is used in " \
           "combination with the assembly mass, the irradiation time and the " \
           "thermal power output to calculate the specific burnup needed for " \
           "the Gaussian process regression." \
  }
  int n_assem_core;

  #pragma cyclus var { \
    "doc": "Number of assemblies constituting a single batch. This is the " \
           "number of assemblies that gets discharged from the core after " \
           "one complete irradiation cycle." \
  }
  int n_assem_batch;

  #pragma cyclus var { \
    "uitype": "range", \
    "range": [1., 1e5], \
    "doc": "The mass of one assembly in kg. This value is used in " \
           "combination with the number of assemblies in a full core, the " \
           "irradiation time and the thermal power output to calculate the " \
           "specific burnup needed for the Gaussian process regression." \
  }
  double assem_size;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Irradiation cycle parameters
  #pragma cyclus var { \
    "default": 0, \
    "doc": "If true, the archetype transmutes all assemblies upon " \
           "decommissioning. If false, the archetype only transmutes half." \
  }
  bool decom_transmute_all;

  #pragma cyclus var { \
    "default": 12, \
    "doc": "The duration of one complete irradiation cycle excluding the " \
           "refuelling, in units of simulation time steps. This value is " \
           "used in combination with the core mass and the thermal power " \
           "output to calculate the specific burnup needed for the Gaussian " \
           "process regression." \
  }
  int cycle_time;

  #pragma cyclus var { \
    "default": 1, \
    "doc": "The duration of an entire refuelling period, i.e., the minimum " \
           "time between the end of a cycle and the start of the next cycle." \
           "In units of simulation time steps." \
  }
  int refuel_time;

  #pragma cyclus var { \
    "default": 0, \
    "doc": "The number of time steps since the start of the last cycle. Only " \
           "set this variable if you know what you are doing." \
  }
  int cycle_step;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Side product settings
  #pragma cyclus var { \
    "default": 1, \
    "internal": True, \
    "doc": "If the reactor produces side products, then this variable is set " \
           "to true." \
  }
  bool is_hybrid;

  #pragma cyclus var { \
    "default": [], \
    "doc": "An ordered vector containing the side products produced during " \
           "normal reactor operation." \
  }
  std::vector<std::string> side_products;

  #pragma cyclus var { \
    "default": [], \
    "doc": "An ordered vector containing the quantities of the produced side " \
           "products. The vector entries must be of type 'double'." \
  }
  std::vector<double> side_product_quantity;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Other variables
  #pragma cyclus var { \
    "default": 0, \
    "uitype": "range", \
    "range": [0, 10000], \
    "units": "MWth", \
    "doc": "The amount of thermal power that the reactor produces during " \
           "operation. This value is used in combination with the core mass " \
           "and the irradiation time to calculate the specific burnup needed " \
           "for the Gaussian process regression." \
  }
  double power_output;

  #pragma cyclus var { \
    "default": 350, \
    "units": "K", \
    "doc": "The reactor moderator temperature." \
  }
  double temperature;

  // This variable is internal only and true if fuel has already been discharged
  // this cycle.
  #pragma cyclus var { \
    "default": 0, \
    "internal": True, \
    "doc": "This variable should NEVER be set manually." \
  }
  bool discharged;

  #pragma cyclus var { \
    "default": 0, \
    "uitype": "range", \
    "range": [0, 1000000000], \
    "doc": "Maximum number of fresh fuel assemblies that are stored in " \
           "the facility if available." \
  }
  int n_assem_fresh;

  #pragma cyclus var { \
    "default": 1000000000, \
    "uitype": "range", \
    "range": [0, 1000000000], \
    "doc": "Maximum number of spent fuel assemblies that can be stored in " \
           "the facility before reactor operation stalls." \
  }
  int n_assem_spent;

  // This variable should be hidden/unavailable in ui.  Maps resource object
  // id's to the index for the incommod through which they were received.
  #pragma cyclus var { \
    "default": {}, \
    "internal": True, \
    "doc": "This variable should NEVER be set manually." \
  }
  std::map<int,int> res_indexes;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Material buffers
  # pragma cyclus var {"capacity": "n_assem_fresh * assem_size"}
  cyclus::toolkit::ResBuf<cyclus::Material> fresh_inv;
  # pragma cyclus var {"capacity": "n_assem_core * assem_size"}
  cyclus::toolkit::ResBuf<cyclus::Material> core;
  # pragma cyclus var {"capacity": "n_assem_spent * assem_size"}
  cyclus::toolkit::ResBuf<cyclus::Material> spent_inv;

  // This set contains all nuc ids that may be part of fresh fuel. So far,
  // they are limited to U235 and U238, however this list will expand in the
  // future. Notably, other U isotopes will be included (probably U232 up to
  // U236) and possibly Pu as well.
  const std::set<int> permitted_fresh_fuel_comps;
  // This set contains all nuc ids of isotopes in spent fuel that we are
  // interested in and that the GPRs calculate.
  const std::set<int> relevant_spent_fuel_comps;

  // Filenames of files used to pass arguments and results between the Python
  // file and the C++ Cyclus archetype.
  std::string out_fname;
  std::string in_fname;

  // This variable stores a constant unique identifier (the time since epoch in
  // ns upon instantiation) and is used when storing the JSON files.
  const uint64_t uid_fname;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Coordinates
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

  bool Discharge_();
  bool Retired_();
  cyclus::Composition::Ptr ImportSpentFuelComposition_(double qty);
  std::map<std::string, cyclus::toolkit::MatVec> PeekSpent_();
  std::map<std::string, cyclus::toolkit::MatVec> PopSpent_();
  std::string OutCommod_(cyclus::Material::Ptr m);
  void CompositionToOutFile_(cyclus::Composition::Ptr comp,
                             bool delete_outfile = false);
  void IndexRes_(cyclus::Resource::Ptr m, std::string incommod);
  void Load_();
  void Record_(std::string name, std::string val);
  void RecordPosition_();
  void RecordSideProduct_(bool produce);
  void PushSpent_(std::map<std::string, cyclus::toolkit::MatVec> mats);
  void Transmute_();
  void Transmute_(int n_assem);

  static uint64_t GetUid_();
};

}  // namespace misoenrichment

#endif  // MISOENRICHMENT_SRC_GPR_REACTOR_H_
