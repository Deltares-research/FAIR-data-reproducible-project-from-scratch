from pathlib import Path

import imod
import numpy as np
import xarray as xr

# %%
# Path management
path_discretization = snakemake.input["path_discretization"]
path_drn_pkg = snakemake.input["path_drn_pkg"]
path_riv_pkg = snakemake.input["path_riv_pkg"]
path_recharge = snakemake.input["path_recharge"]
path_ic = snakemake.input["path_ic"]
path_chd = snakemake.input["path_chd"]
path_subsurface = snakemake.input["path_subsurface"]

path_model = snakemake.output["path_model"]
path_model_directory = Path(path_model).parent

# %% Fill model

gwf_model = imod.mf6.GroundwaterFlowModel()
gwf_model["dis"] = imod.mf6.StructuredDiscretization.from_file(path_discretization)
gwf_model["drn"] = imod.mf6.Drainage.from_file(path_drn_pkg)
gwf_model["riv"] = imod.mf6.River.from_file(path_riv_pkg)
gwf_model["rch"] = imod.mf6.Recharge.from_file(path_recharge)
gwf_model["ic"] = imod.mf6.InitialConditions.from_file(path_ic)
gwf_model["chd"] = imod.mf6.ConstantHead.from_file(path_chd)
gwf_model["npf"] = imod.mf6.NodePropertyFlow.from_file(path_subsurface)

ss = 0.0003
layer = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
sy = xr.DataArray(
    [0.16, 0.16, 0.16, 0.16, 0.15, 0.15, 0.15, 0.15, 0.14, 0.14, 0.14, 0.14, 0.14],
    {"layer": layer},
    ("layer",),
)
times_sto = np.array(
    [
        "2009-12-30T23:59:59.00",
        "2009-12-31T00:00:00.00",
        "2010-12-31T00:00:00.00",
        "2011-12-31T00:00:00.00",
        "2012-12-31T00:00:00.00",
        "2013-12-31T00:00:00.00",
        "2014-12-31T00:00:00.00",
    ],
    dtype="datetime64[ns]",
)

transient = xr.DataArray(
    [False, True, True, True, True, True, True], {"time": times_sto}, ("time",)
)

gwf_model["sto"] = imod.mf6.SpecificStorage(
    specific_storage=ss,
    specific_yield=sy,
    transient=transient,
    convertible=0,
    save_flows=True,
)

gwf_model["oc"] = imod.mf6.OutputControl(save_head="last", save_budget="last")

simulation = imod.mf6.Modflow6Simulation("hondsrug_model")
simulation["GWF"] = gwf_model

# %%
# Solver settings
# ---------------
#
# The solver settings are indicated using
# :doc:`/api/generated/mf6/imod.mf6.Solution`.
# If the values are not indicated manually, the defaults values will be considered.

simulation["solver"] = imod.mf6.Solution(
    modelnames=["GWF"],
    print_option="summary",
    csv_output=False,
    no_ptc=True,
    outer_dvclose=1.0e-4,
    outer_maximum=500,
    under_relaxation=None,
    inner_dvclose=1.0e-4,
    inner_rclose=0.001,
    inner_maximum=100,
    linear_acceleration="cg",
    scaling_method=None,
    reordering_method=None,
    relaxation_factor=0.97,
)

# %%
# Assign time discretization
# --------------------------
#
# The time discretization of this model is 6 years.

simulation.create_time_discretization(
    additional_times=["2009-12-30T23:59:59.000000000", "2015-12-31T00:00:00.000000000"]
)

simulation.write(path_model_directory, binary=False)