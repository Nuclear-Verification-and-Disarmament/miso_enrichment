import sys

sys.path.append("./")
sys.path.append("./input")

import archetypes
import commodity
import control
import facility
import recipe
import region


def simulation():
    arch = archetypes.archetypes()
    commod = commodity.commodity()
    ctrl = control.control()
    fac = facility.facility()
    recipes = recipe.recipe()
    reg = region.region()

    return {"simulation": {**arch, **commod, **ctrl, **fac, **recipes, **reg}}
