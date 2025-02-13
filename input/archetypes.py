def archetypes():
    d = {
        "archetypes": {
            "spec": [
                {"lib": "agents", "name": "NullInst"},
                {"lib": "agents", "name": "NullRegion"},
                {"lib": "cycamore", "name": "Sink"},
                {"lib": "cycamore", "name": "Source"},
                {"lib": "cycamore", "name": "Storage"},
                #{"lib": "misoenrichment", "name": "GprReactor"},
                {"lib": "misoenrichment", "name": "MIsoEnrich"},
                {"lib": "misoenrichment", "name": "VarRecipeSource"},
            ]
        }
    }

    return d
