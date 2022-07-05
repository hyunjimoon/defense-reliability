"""
Python model 'fixedQuantity.py'
Translated using PySD
"""

from pathlib import Path
import numpy as np
import xarray as xr
from scipy import stats

from pysd.py_backend.functions import if_then_else
from pysd.py_backend.statefuls import Smooth, Integ
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
    name="FINAL TIME", units="Day", comp_type="Constant", comp_subtype="Normal"
)
def final_time():
    """
    The final time for the simulation.
    """
    return __data["time"].final_time()


@component.add(
    name="INITIAL TIME", units="Day", comp_type="Constant", comp_subtype="Normal"
)
def initial_time():
    """
    The initial time for the simulation.
    """
    return __data["time"].initial_time()


@component.add(
    name="SAVEPER",
    units="Day",
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
    units="Day",
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


@component.add(name="Fixed Order Quantity", comp_type="Constant", comp_subtype="Normal")
def fixed_order_quantity():
    return 1000


@component.add(
    name="Supply Line",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_supply_line": 1},
    other_deps={
        "_integ_supply_line": {
            "initial": {"demand_forecast": 1, "perceived_lead_time": 1},
            "step": {"purchasing": 1, "shipping": 1},
        }
    },
)
def supply_line():
    return _integ_supply_line()


_integ_supply_line = Integ(
    lambda: purchasing() - shipping(),
    lambda: demand_forecast() * perceived_lead_time(),
    "_integ_supply_line",
)


@component.add(
    name="Service Policy Coefficient", comp_type="Constant", comp_subtype="Normal"
)
def service_policy_coefficient():
    return 2


@component.add(
    name="Purchase Amount",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "supply_line": 1,
        "inventory": 1,
        "reorder_quantity": 1,
        "fixed_order_quantity": 1,
    },
)
def purchase_amount():
    return if_then_else(
        supply_line() + inventory() <= reorder_quantity(),
        lambda: fixed_order_quantity(),
        lambda: 0,
    )


@component.add(
    name="Reorder Quantity",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "demand_forecast": 1,
        "perceived_lead_time": 2,
        "sd_demand": 1,
        "service_policy_coefficient": 1,
    },
)
def reorder_quantity():
    return (
        demand_forecast() * perceived_lead_time()
        + service_policy_coefficient() * sd_demand() * np.sqrt(perceived_lead_time())
    )


@component.add(
    name="Inventory",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_inventory": 1},
    other_deps={
        "_integ_inventory": {
            "initial": {"fixed_order_quantity": 1},
            "step": {"shipping": 1, "selling": 1},
        }
    },
)
def inventory():
    return _integ_inventory()


_integ_inventory = Integ(
    lambda: shipping() - selling(), lambda: fixed_order_quantity(), "_integ_inventory"
)


@component.add(name="Lead Time Adjustment", comp_type="Constant", comp_subtype="Normal")
def lead_time_adjustment():
    return 5


@component.add(
    name="Underage Cost",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"deficiency": 1, "unit_underage_cost": 1},
)
def underage_cost():
    return deficiency() * unit_underage_cost()


@component.add(
    name="Cost",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"underage_cost": 1, "overage_cost": 1},
)
def cost():
    return underage_cost() + overage_cost()


@component.add(name="Unit Underage Cost", comp_type="Constant", comp_subtype="Normal")
def unit_underage_cost():
    return 9


@component.add(name="Unit Overage Cost", comp_type="Constant", comp_subtype="Normal")
def unit_overage_cost():
    return 1


@component.add(name="Mean Demand", comp_type="Constant", comp_subtype="Normal")
def mean_demand():
    return 100


@component.add(name="Mean Lead Time", comp_type="Constant", comp_subtype="Normal")
def mean_lead_time():
    return 5


@component.add(name="Std Lead Time", comp_type="Constant", comp_subtype="Normal")
def std_lead_time():
    return 3


@component.add(
    name="Overage Cost",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"inventory": 1, "supply_line": 1, "unit_overage_cost": 1},
)
def overage_cost():
    return (inventory() + supply_line()) * unit_overage_cost()


@component.add(name="Sd Demand", comp_type="Constant", comp_subtype="Normal")
def sd_demand():
    return 50


@component.add(
    name="Perceived Lead Time",
    comp_type="Stateful",
    comp_subtype="Smooth",
    depends_on={"_smooth_perceived_lead_time": 1},
    other_deps={
        "_smooth_perceived_lead_time": {
            "initial": {"lead_time": 1, "lead_time_adjustment": 1},
            "step": {"lead_time": 1, "lead_time_adjustment": 1},
        }
    },
)
def perceived_lead_time():
    return _smooth_perceived_lead_time()


_smooth_perceived_lead_time = Smooth(
    lambda: lead_time(),
    lambda: lead_time_adjustment(),
    lambda: lead_time(),
    lambda: 1,
    "_smooth_perceived_lead_time",
)


@component.add(
    name="Deficiency",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"backlog": 1, "selling": 1},
)
def deficiency():
    return np.maximum(0, backlog() - selling())


@component.add(
    name="Backlog Out",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"selling": 1},
)
def backlog_out():
    return selling()


@component.add(
    name="Backlog",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_backlog": 1},
    other_deps={
        "_integ_backlog": {"initial": {}, "step": {"backlog_in": 1, "backlog_out": 1}}
    },
)
def backlog():
    return _integ_backlog()


_integ_backlog = Integ(
    lambda: backlog_in() - backlog_out(), lambda: 100, "_integ_backlog"
)


@component.add(
    name="Selling",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"backlog": 1, "inventory": 1},
)
def selling():
    return np.minimum(backlog(), inventory())


@component.add(
    name="Backlog In",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"demand": 1},
)
def backlog_in():
    return demand()


@component.add(
    name="Purchasing",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"purchase_amount": 1},
)
def purchasing():
    return purchase_amount()


@component.add(
    name="Lead Time",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"mean_lead_time": 1, "std_lead_time": 1},
)
def lead_time():
    return stats.truncnorm.rvs(
        1, 10, loc=mean_lead_time(), scale=std_lead_time(), size=()
    )


@component.add(
    name="Shipping",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"supply_line": 1, "lead_time": 1},
)
def shipping():
    return supply_line() / lead_time()


@component.add(
    name="Demand",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"mean_demand": 1, "sd_demand": 1},
)
def demand():
    return stats.truncnorm.rvs(0, 200, loc=mean_demand(), scale=sd_demand(), size=())


@component.add(
    name="Demand Forecast",
    comp_type="Stateful",
    comp_subtype="Smooth",
    depends_on={"_smooth_demand_forecast": 1},
    other_deps={
        "_smooth_demand_forecast": {
            "initial": {"demand": 1, "forecaset_period": 1},
            "step": {"demand": 1, "forecaset_period": 1},
        }
    },
)
def demand_forecast():
    return _smooth_demand_forecast()


_smooth_demand_forecast = Smooth(
    lambda: demand(),
    lambda: forecaset_period(),
    lambda: demand(),
    lambda: 1,
    "_smooth_demand_forecast",
)


@component.add(name="Forecaset Period", comp_type="Constant", comp_subtype="Normal")
def forecaset_period():
    return 5
