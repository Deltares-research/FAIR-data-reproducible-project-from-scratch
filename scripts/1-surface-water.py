import imod
import xarray as xr

# %%
# Path management
path_drainage = snakemake.input["path_drainage"]
path_river = snakemake.input["path_river"]
path_drn_pkg = snakemake.output["path_drn_pkg"]
path_riv_pkg = snakemake.output["path_riv_pkg"]

# %%
# Open dataset
drainage = xr.open_dataset(path_drainage)
river = xr.open_dataset(path_river)

# %%
# Drainage package - DRN
# =======================
#
# The drain package is used to simulate features that remove water from the aquifer,
# such as agricultural drains or springs.
# This occurs at a rate proportional to the head difference between the head in the
# aquifer and the drain elevation
# (the head in the aquifer has to be above that elevation).
# The conductance is the proportionality constant.
#
# Import drainage information
# ----------------------------

pipe_cond = drainage["conductance"]
pipe_elev = drainage["elevation"]

pipe_cond

# %%
# Adding information to the DRN package
# -------------------------------------
#
# To add the information to the DRN package within the gwf_model variable, the
# :doc:`/api/generated/mf6/imod.mf6.Drainage`.
# function is used. It is required to add the previously created arrays for
# the drain elevation and the drain conductance.
# It is optional to insert the information for
# ``print_input``, ``print_flows`` and ``save_flows``
# which are set to False by default.

drn_pkg = imod.mf6.Drainage(conductance=pipe_cond, elevation=pipe_elev)

# %%
# River package - RIV
# ===================
#
# This package simulates the effects of flow between
# surface-water features and groundwater systems.
#
# Import river information
# ------------------------

riv_cond = river["conductance"]
riv_stage = river["stage"]
riv_bot = river["bottom"]

# %%
# Adding information to the RIV package
# -------------------------------------
#
# The data is assigned to the gwf_model variable by using
# :doc:`/api/generated/mf6/imod.mf6.River`,
# based on the previously imported conductance, stage and bottom arrays.

riv_pkg = imod.mf6.River(
    conductance=riv_cond, stage=riv_stage, bottom_elevation=riv_bot
)

# %%
riv_pkg.to_netcdf(path_riv_pkg)
drn_pkg.to_netcdf(path_drn_pkg)