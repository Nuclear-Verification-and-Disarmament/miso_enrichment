import facility


def institution():
    fac = facility.facility()
    d = {
        "institution": [
            {
                "name": "MyInstitution",
                "config": {"NullInst": None},
                "initialfacilitylist": {
                    "entry": [
                        {"number": 1, "prototype": f["name"]} for f in fac["facility"]
                    ]
                },
            }
        ]
    }
    return d
