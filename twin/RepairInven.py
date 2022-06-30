"""
Python model 'RepairInven.py'
Translated using PySD
"""

from pathlib import Path
import numpy as np
import xarray as xr

from pysd.py_backend.functions import if_then_else, modulo, xidz
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


@component.add(name="Purchasing Period", comp_type="Constant", comp_subtype="Normal")
def purchasing_period():
    return 12


@component.add(
    name="BackLog",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_backlog": 1},
    other_deps={
        "_integ_backlog": {"initial": {}, "step": {"backlogin": 1, "backlogout": 1}}
    },
)
def backlog():
    return _integ_backlog()


_integ_backlog = Integ(lambda: backlogin() - backlogout(), lambda: 0, "_integ_backlog")


@component.add(
    name="BackLogIn",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"demand": 1},
)
def backlogin():
    return demand()


@component.add(
    name="BackLogOut",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"using": 1},
)
def backlogout():
    return using()


@component.add(
    name="Battle Field",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_battle_field": 1},
    other_deps={
        "_integ_battle_field": {
            "initial": {"initial_value": 1},
            "step": {"input_1": 1, "maintenance": 1},
        }
    },
)
def battle_field():
    return _integ_battle_field()


_integ_battle_field = Integ(
    lambda: input_1() - maintenance(), lambda: initial_value(), "_integ_battle_field"
)


@component.add(
    name="Cost",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"deficiency_cost": 1, "inventory_cost": 1},
)
def cost():
    return deficiency_cost() + inventory_cost()


@component.add(
    name="Deficiency",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"backlog": 1, "maintenance_time": 1, "using": 1},
)
def deficiency():
    return np.maximum(0, backlog() / maintenance_time() - using())


@component.add(
    name="Deficiency Cost",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"deficiency_unit_cost": 1, "deficiency": 1},
)
def deficiency_cost():
    return deficiency_unit_cost() * deficiency()


@component.add(name="Deficiency Unit Cost", comp_type="Constant", comp_subtype="Normal")
def deficiency_unit_cost():
    return 99


@component.add(
    name="Demand",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"maintenance": 1},
)
def demand():
    return maintenance()


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


@component.add(name="Failure Rate", comp_type="Constant", comp_subtype="Normal")
def failure_rate():
    return 0.01


@component.add(name="Forecast Period", comp_type="Constant", comp_subtype="Normal")
def forecast_period():
    return 3


@component.add(name="Initial Value", comp_type="Constant", comp_subtype="Normal")
def initial_value():
    return 100


@component.add(
    name="Input",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"repair_shop": 1, "maintenance_time": 1, "inventory": 1, "backlog": 1},
)
def input_1():
    return (
        repair_shop()
        / maintenance_time()
        * np.minimum(1, xidz(inventory(), backlog(), 1))
    )


@component.add(
    name="Inventory",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_inventory": 1},
    other_deps={
        "_integ_inventory": {
            "initial": {"target_inventory": 1},
            "step": {"shipping": 1, "using": 1},
        }
    },
)
def inventory():
    return _integ_inventory()


_integ_inventory = Integ(
    lambda: shipping() - using(), lambda: target_inventory(), "_integ_inventory"
)


@component.add(
    name="Inventory Adjustment Value",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"target_inventory": 1, "inventory": 1},
)
def inventory_adjustment_value():
    return target_inventory() - inventory()


@component.add(
    name="Inventory Cost",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"inventory": 1, "supply_line": 1, "inventory_unit_cost": 1},
)
def inventory_cost():
    return (inventory() + supply_line()) * inventory_unit_cost()


@component.add(name="Inventory Unit Cost", comp_type="Constant", comp_subtype="Normal")
def inventory_unit_cost():
    return 1


@component.add(name="Lead Time", comp_type="Constant", comp_subtype="Normal")
def lead_time():
    return 10


@component.add(
    name="Purchasing",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"purchase_amount": 1},
)
def purchasing():
    return purchase_amount()


@component.add(name="Maintenance Time", comp_type="Constant", comp_subtype="Normal")
def maintenance_time():
    return 1


@component.add(
    name="Shipping",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"supply_line": 1, "lead_time": 1},
)
def shipping():
    return supply_line() / lead_time()


@component.add(name="N of Days in Stock", comp_type="Constant", comp_subtype="Normal")
def n_of_days_in_stock():
    return 2


@component.add(
    name="Supply Line",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_supply_line": 1},
    other_deps={
        "_integ_supply_line": {
            "initial": {"target_supply_line": 1},
            "step": {"purchasing": 1, "shipping": 1},
        }
    },
)
def supply_line():
    return _integ_supply_line()


_integ_supply_line = Integ(
    lambda: purchasing() - shipping(),
    lambda: target_supply_line(),
    "_integ_supply_line",
)


@component.add(
    name="Predictive Maintenance Period", comp_type="Constant", comp_subtype="Normal"
)
def predictive_maintenance_period():
    return 5


@component.add(
    name="Purchase Adjustment Value",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"target_supply_line": 1, "supply_line": 1},
)
def purchase_adjustment_value():
    return target_supply_line() - supply_line()


@component.add(
    name="Purchase Amount",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"time": 1, "purchasing_period": 1, "inventory_adjustment_value": 1},
)
def purchase_amount():
    return if_then_else(
        modulo(time(), purchasing_period()) == 0,
        lambda: inventory_adjustment_value(),
        lambda: 0,
    )


@component.add(
    name="Target Supply Line",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"lead_time": 1, "demand_forecast": 1},
)
def target_supply_line():
    return lead_time() * demand_forecast()


@component.add(
    name="Repair Shop",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_repair_shop": 1},
    other_deps={
        "_integ_repair_shop": {"initial": {}, "step": {"maintenance": 1, "input_1": 1}}
    },
)
def repair_shop():
    return _integ_repair_shop()


_integ_repair_shop = Integ(
    lambda: maintenance() - input_1(), lambda: 0, "_integ_repair_shop"
)


@component.add(
    name="Target Inventory",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"demand_forecast": 1, "n_of_days_in_stock": 1},
)
def target_inventory():
    return demand_forecast() * n_of_days_in_stock()


@component.add(
    name="Using",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"backlog": 1, "input_1": 1},
)
def using():
    return np.minimum(backlog(), input_1())


@component.add(
    name="Maintenance",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"n_of_contigency_maintenance": 1, "n_of_predictive_maintenance": 1},
)
def maintenance():
    return n_of_contigency_maintenance() + n_of_predictive_maintenance()


@component.add(
    name="N of Contigency Maintenance",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"battle_field": 1, "n_of_predictive_maintenance": 1, "failure_rate": 1},
)
def n_of_contigency_maintenance():
    return (battle_field() - n_of_predictive_maintenance()) * failure_rate()


@component.add(
    name="N of Predictive Maintenance",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"battle_field": 1, "predictive_maintenance_period": 1},
)
def n_of_predictive_maintenance():
    return battle_field() / predictive_maintenance_period()
