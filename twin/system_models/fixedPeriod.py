"""
Python model 'fixedPeriod.py'
Translated using PySD
"""

from pathlib import Path
import numpy as np
import xarray as xr
from scipy import stats

from pysd.py_backend.functions import modulo, if_then_else
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


@component.add(name="Fixed Period", comp_type="Constant", comp_subtype="Normal")
def fixed_period():
    return 7


@component.add(
    name="Inventory",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_inventory": 1},
    other_deps={
        "_integ_inventory": {
            "initial": {"target_inventory": 1, "mean_lead_time": 1},
            "step": {"shipping": 1, "sales": 1},
        }
    },
)
def inventory():
    return _integ_inventory()


_integ_inventory = Integ(
    lambda: shipping() - sales(),
    lambda: target_inventory() - mean_lead_time() * 100,
    "_integ_inventory",
)


@component.add(
    name="Purchase Amount",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "time": 1,
        "fixed_period": 2,
        "target_inventory": 1,
        "inventory": 1,
        "supplyline": 1,
    },
)
def purchase_amount():
    return if_then_else(
        modulo(time(), fixed_period()) == fixed_period() - 1,
        lambda: np.maximum(0, target_inventory() - supplyline() - inventory()),
        lambda: 0,
    )


@component.add(name="Target Inventory", comp_type="Constant", comp_subtype="Normal")
def target_inventory():
    return 1500


@component.add(
    name="SupplyLine",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_supplyline": 1},
    other_deps={
        "_integ_supplyline": {
            "initial": {"demand_forecast": 1, "perceived_lead_time": 1},
            "step": {"purchasing": 1, "shipping": 1},
        }
    },
)
def supplyline():
    return _integ_supplyline()


_integ_supplyline = Integ(
    lambda: purchasing() - shipping(),
    lambda: demand_forecast() * perceived_lead_time(),
    "_integ_supplyline",
)


@component.add(
    name="Lead Time Adjustment Period", comp_type="Constant", comp_subtype="Normal"
)
def lead_time_adjustment_period():
    return 5


@component.add(
    name="Underage Cost",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"deficiency": 1, "underage_unit_cost": 1},
)
def underage_cost():
    return deficiency() * underage_unit_cost()


@component.add(
    name="Cost",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"underage_cost": 1, "overage_cost": 1},
)
def cost():
    return underage_cost() + overage_cost()


@component.add(name="Underage Unit Cost", comp_type="Constant", comp_subtype="Normal")
def underage_unit_cost():
    return 9


@component.add(name="Overage Unit Cost", comp_type="Constant", comp_subtype="Normal")
def overage_unit_cost():
    return 1


@component.add(name="Mean Demand", comp_type="Constant", comp_subtype="Normal")
def mean_demand():
    return 100


@component.add(name="Mean Lead Time", comp_type="Constant", comp_subtype="Normal")
def mean_lead_time():
    return 5


@component.add(name="Sd Lead Time", comp_type="Constant", comp_subtype="Normal")
def sd_lead_time():
    return 3


@component.add(
    name="Overage Cost",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"inventory": 1, "supplyline": 1, "overage_unit_cost": 1},
)
def overage_cost():
    return (inventory() + supplyline()) * overage_unit_cost()


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
            "initial": {"lead_time": 1, "lead_time_adjustment_period": 1},
            "step": {"lead_time": 1, "lead_time_adjustment_period": 1},
        }
    },
)
def perceived_lead_time():
    return _smooth_perceived_lead_time()


_smooth_perceived_lead_time = Smooth(
    lambda: lead_time(),
    lambda: lead_time_adjustment_period(),
    lambda: lead_time(),
    lambda: 1,
    "_smooth_perceived_lead_time",
)


@component.add(
    name="Deficiency",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"backlog": 1, "sales": 1},
)
def deficiency():
    return np.maximum(0, backlog() - sales())


@component.add(
    name="BacklogOut",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"sales": 1},
)
def backlogout():
    return sales()


@component.add(
    name="Backlog",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_backlog": 1},
    other_deps={
        "_integ_backlog": {"initial": {}, "step": {"backlogin": 1, "backlogout": 1}}
    },
)
def backlog():
    return _integ_backlog()


_integ_backlog = Integ(
    lambda: backlogin() - backlogout(), lambda: 100, "_integ_backlog"
)


@component.add(
    name="Sales",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"backlog": 1, "inventory": 1},
)
def sales():
    return np.minimum(backlog(), inventory())


@component.add(
    name="BacklogIn",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"demand": 1},
)
def backlogin():
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
    depends_on={"mean_lead_time": 1, "sd_lead_time": 1},
)
def lead_time():
    return stats.truncnorm.rvs(
        1, 10, loc=mean_lead_time(), scale=sd_lead_time(), size=()
    )


@component.add(
    name="Shipping",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"supplyline": 1, "lead_time": 1},
)
def shipping():
    return supplyline() / lead_time()


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


@component.add(name="Forecast Period", comp_type="Constant", comp_subtype="Normal")
def forecast_period():
    return 5
