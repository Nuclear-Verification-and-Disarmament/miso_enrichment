def recipe():
    d = {
        "recipe": [
            {
                "name": "NaturalURecipe",
                "basis": "mass",
                "nuclide": [
                    # At the moment (May 2021), U234 is not taken into
                    # account by the Gpr model.
                    # {"id": "U234", "comp": 5.5e-3},
                    {"id": "U235", "comp": 0.711},
                    {"id": "U238", "comp": 100 - 0.711},
                ],
            },
            {
                "name": "EnrichedURecipe",
                "basis": "mass",
                "nuclide": [{"id": "U235", "comp": 1.1}, {"id": "U238", "comp": 98.9}],
            },
            {
                "name": "DepletedURecipe",
                "basis": "mass",
                "nuclide": [{"id": "U235", "comp": 0.3}, {"id": "U238", "comp": 99.7}],
            },
        ]
    }
    return d
