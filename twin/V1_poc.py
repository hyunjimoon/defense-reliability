"""
Python model 'V1_poc.py'
Translated using PySD
"""

from pathlib import Path
import numpy as np
import xarray as xr
from scipy import stats

from pysd.py_backend.statefuls import Integ, Smooth
from pysd import Component

__pysd_version__ = "3.3.0"

__data = {"scope": None, "time": lambda: 0}

_root = Path(__file__).parent


component = Component()

#######################################################################
#                          CONTROL VARIABLES                          #
#######################################################################

_control_vars = {
    "initial_time": lambda: 0,
    "final_time": lambda: 100,
    "time_step": lambda: 1,
    "saveper": lambda: time_step(),
}


def _init_outer_references(data):
    for key in data:
        __data[key] = data[key]


@component.add(name="Time")
def time():
    """
    Current time of the model.
    """
    return __data["time"]()


@component.add(
    name="FINAL TIME", units="Month", comp_type="Constant", comp_subtype="Normal"
)
def final_time():
    """
    The final time for the simulation.
    """
    return __data["time"].final_time()


@component.add(
    name="INITIAL TIME", units="Month", comp_type="Constant", comp_subtype="Normal"
)
def initial_time():
    """
    The initial time for the simulation.
    """
    return __data["time"].initial_time()


@component.add(
    name="SAVEPER",
    units="Month",
    limits=(0.0, np.nan),
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"time_step": 1},
)
def saveper():
    """
    The frequency with which output is stored.
    """
    return __data["time"].saveper()


@component.add(
    name="TIME STEP",
    units="Month",
    limits=(0.0, np.nan),
    comp_type="Constant",
    comp_subtype="Normal",
)
def time_step():
    """
    The time step for the simulation.
    """
    return __data["time"].time_step()


#######################################################################
#                           MODEL VARIABLES                           #
#######################################################################


@component.add(
    name="Adjustment for Inventory",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"desired_inventory": 1, "inventory": 1},
)
def adjustment_for_inventory():
    return desired_inventory() - inventory()


@component.add(
    name="Adjustment for Supply Line",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"desired_supply_line": 1, "supply_line": 1},
)
def adjustment_for_supply_line():
    return desired_supply_line() - supply_line()


@component.add(
    name="Back Log",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_back_log": 1},
    other_deps={"_integ_back_log": {"initial": {}, "step": {"bl_in": 1, "bl_out": 1}}},
)
def back_log():
    return _integ_back_log()


_integ_back_log = Integ(lambda: bl_in() - bl_out(), lambda: 100, "_integ_back_log")


@component.add(
    name="BL In", comp_type="Auxiliary", comp_subtype="Normal", depends_on={"demand": 1}
)
def bl_in():
    return demand()


@component.add(
    name="BL Out",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"shipment_rate": 1},
)
def bl_out():
    return shipment_rate()


@component.add(
    name="Cost",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"underage_cost": 1, "overage_cost": 1},
)
def cost():
    return underage_cost() + overage_cost()


@component.add(
    name="Deficient Amount",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"back_log": 1, "shipment_rate": 1},
)
def deficient_amount():
    return np.maximum(0, back_log() - shipment_rate())


@component.add(
    name="Delivery Rate",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"supply_line": 1, "lead_time": 1},
)
def delivery_rate():
    return supply_line() / lead_time()


@component.add(
    name="Demand",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"mean_of_demand": 1, "sd_of_demand": 1},
)
def demand():
    return stats.truncnorm.rvs(
        0, 200, loc=mean_of_demand(), scale=sd_of_demand(), size=()
    )


@component.add(
    name="Demand Forecast",
    comp_type="Stateful",
    comp_subtype="Smooth",
    depends_on={"_smooth_demand_forecast": 1},
    other_deps={
        "_smooth_demand_forecast": {
            "initial": {"demand": 1, "forecast_period": 1},
            "step": {"demand": 1, "forecast_period": 1},
        }
    },
)
def demand_forecast():
    return _smooth_demand_forecast()


_smooth_demand_forecast = Smooth(
    lambda: demand(),
    lambda: forecast_period(),
    lambda: demand(),
    lambda: 1,
    "_smooth_demand_forecast",
)


@component.add(
    name="Desired Inventory",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"demand_forecast": 1, "inventory_period": 1},
)
def desired_inventory():
    return demand_forecast() * inventory_period()


@component.add(
    name="Desired Supply Line",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"demand_forecast": 1, "lead_time": 1},
)
def desired_supply_line():
    return demand_forecast() * lead_time()


@component.add(name="Forecast Period", comp_type="Constant", comp_subtype="Normal")
def forecast_period():
    return 3


@component.add(
    name="Inventory",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_inventory": 1},
    other_deps={
        "_integ_inventory": {
            "initial": {"desired_inventory": 1},
            "step": {"delivery_rate": 1, "shipment_rate": 1},
        }
    },
)
def inventory():
    return _integ_inventory()


_integ_inventory = Integ(
    lambda: delivery_rate() - shipment_rate(),
    lambda: desired_inventory(),
    "_integ_inventory",
)


@component.add(name="Inventory Period", comp_type="Constant", comp_subtype="Normal")
def inventory_period():
    return 5


@component.add(name="Lead Time", comp_type="Constant", comp_subtype="Normal")
def lead_time():
    return 5


@component.add(name="Mean of Demand", comp_type="Constant", comp_subtype="Normal")
def mean_of_demand():
    return 100


@component.add(
    name="Overage Cost",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"inventory": 1, "supply_line": 1, "unit_overage_cost": 1},
)
def overage_cost():
    return (inventory() + supply_line()) * unit_overage_cost()


@component.add(
    name="Purchase",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "adjustment_for_inventory": 1,
        "adjustment_for_supply_line": 1,
        "demand_forecast": 1,
    },
)
def purchase():
    return np.maximum(
        0, adjustment_for_inventory() + adjustment_for_supply_line() + demand_forecast()
    )


@component.add(
    name="Purchase Rate",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"purchase": 1},
)
def purchase_rate():
    return purchase()


@component.add(name="Sd of Demand", comp_type="Constant", comp_subtype="Normal")
def sd_of_demand():
    return 50


@component.add(
    name="Shipment Rate",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"back_log": 1, "inventory": 1},
)
def shipment_rate():
    return np.minimum(back_log(), inventory())


@component.add(
    name="Supply Line",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_supply_line": 1},
    other_deps={
        "_integ_supply_line": {
            "initial": {"desired_supply_line": 1},
            "step": {"purchase_rate": 1, "delivery_rate": 1},
        }
    },
)
def supply_line():
    return _integ_supply_line()


_integ_supply_line = Integ(
    lambda: purchase_rate() - delivery_rate(),
    lambda: desired_supply_line(),
    "_integ_supply_line",
)


@component.add(
    name="Underage Cost",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"deficient_amount": 1, "unit_underage_cost": 1},
)
def underage_cost():
    return deficient_amount() * unit_underage_cost()


@component.add(name="Unit Overage Cost", comp_type="Constant", comp_subtype="Normal")
def unit_overage_cost():
    return 1


@component.add(name="Unit Underage Cost", comp_type="Constant", comp_subtype="Normal")
def unit_underage_cost():
    return 9
