def archetypes():
    d = {"archetypes": {
           "spec": [
             {"lib": "agents", "name": "NullInst"},
             {"lib": "agents", "name": "NullRegion"},
             {"lib": "cycamore", "name": "Sink"},
             {"lib": "cycamore", "name": "Source"},
             {"lib": "misoenrichment", "name":"MIsoEnrich"}
           ]
         }}

    return d
