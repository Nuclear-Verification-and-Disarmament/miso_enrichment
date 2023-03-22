import institution


def region():
    inst = institution.institution()
    d = {
        "region": [
            {
                "name": "MyRegion",
                "config": {"NullRegion": None},
                "institution": inst["institution"],
            }
        ]
    }
    return d
