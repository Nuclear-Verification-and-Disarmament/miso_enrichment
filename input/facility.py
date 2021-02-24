def facility():
    d = {"facility": [
          {
            "name": "NaturalUSource",
            "config": {"Source": {
              "outcommod": "NaturalU",
              "outrecipe": "NaturalURecipe",
              "throughput": 2000
            }}
          },
          {
            "name": "EnrichedUSink",
            "config": {"Sink": {
              "in_commods": {"val": ["EnrichedU"]},
              "recipe_name": "EnrichedURecipe"
            }}
          },
          {
            "name": "DepletedUSink",
            "config": {"Sink": {
              "in_commods": {"val": ["DepletedU"]},
              "recipe_name": "DepletedURecipe"
            }}
          },
          {
            "name": "EnrichmentFacility",
            "config": {"MIsoEnrich": {
              "feed_commod": "NaturalU",
              "feed_recipe": "NaturalURecipe",
              "product_commod": "EnrichedU",
              "tails_commod": "DepletedU",
              "tails_assay": 0.003,
              "initial_feed": 0,
              "max_feed_inventory": 10000,
              "gamma_235": 1.35,
              "swu_capacity": 1e299,
              "swu_capacity_vals": {"val": [1e5, 1000, 2000]},
              "swu_capacity_times": {"val": [0, 5, 6]},
              "use_downblending": True
            }}
          }
        ]}
    return d
