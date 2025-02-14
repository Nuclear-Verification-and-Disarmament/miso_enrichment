def facility():
    d = {
        "facility": [
            {
                "name": "NaturalUSource",
                "config": {
                    "VarRecipeSource": {
                        "out_commod": "NaturalU",
                        "var_out_recipe": {
                            "mass_or_atom": "mass",
                            "nuclides": {"item": [
                                {
                                    "nuc_id": 922340000,
                                    "rng_properties": {
                                        "distribution": "uniform",
                                        "parameters": {"val": [5.1e-5, 5.4e-5]}
                                    }
                                },
                                {
                                    "nuc_id": 922350000,
                                    "rng_properties": {
                                        "distribution": "normal",
                                        "parameters": {"val": [7.1e-3, 0.01e-3, 0., 1.]}
                                    }
                                },
                                {
                                    "nuc_id": 922380000,
                                    "rng_properties": {
                                        "distribution": "normalisation",
                                        "parameters": {"val": [0., 0.]}
                                    }
                                }
                            ]}
                        },
                        "throughput_vals": {"val": [10]},
                        "throughput_times": {"val": [0]}
                    }
                }
            },
            {
                "name": "SpentFuelSink",
                "config": {"Sink": {"in_commods": {"val": ["SpentFuel"]}}},
            },
            {
                "name": "DepletedUSink",
                "config": {
                    "Sink": {
                        "in_commods": {"val": ["DepletedU"]},
                        "recipe_name": "DepletedURecipe",
                    }
                },
            },
            {
                "name": "FreshFuelStorage",
                "config": {
                    "Storage": {
                        "in_commods": {"val": ["EnrichedU"]},
                        "out_commods": {"val": ["FreshFuel"]},
                        "in_recipe": ["EnrichedURecipe"],
                        "residence_time": 0,
                    }
                },
            },
            {
                "name": "EnrichmentFacility",
                "config": {
                    "MIsoEnrich": {
                        "feed_commod": "NaturalU",
                        "feed_recipe": "NaturalURecipe",
                        "product_commod": "EnrichedU",
                        "tails_commod": "DepletedU",
                        "tails_assay": 0.003,
                        "initial_feed": 0,
                        "max_feed_inventory": 1e299,
                        "gamma_235": 1.35,
                        "swu_capacity": 1e299,
                        "swu_capacity_vals": {"val": [1e5, 5e4, 5e5]},
                        "swu_capacity_times": {"val": [0, 5, 6]},
                        "use_downblending": True,
                        "use_integer_stages": True,  # Must be 'True' because downblending is enabled
                    }
                },
            },
            #{
            #    "name": "SavannahRiverReactor",
            #    "config": {
            #        "GprReactor": {
            #            "in_commods": {"val": ["FreshFuel"]},
            #            "out_commods": {"val": ["SpentFuel"]},
            #            "in_recipes": {"val": ["EnrichedURecipe"]},
            #            "n_assem_core": 1,
            #            "n_assem_batch": 1,
            #            "assem_size": 110820,
            #            "cycle_time": 88,
            #            "refuel_time": 6,
            #            "power_output": 2400,
            #            "temperature": 350,
            #        }
            #    },
            #},
        ]
    }
    return d
