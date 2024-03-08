import imod
import xarray as xr

# Path management
path_layermodel = snakemake.input["path_layermodel"]
path_subsurface = snakemake.output["path_subsurface"]

# %%
layermodel = xr.open_dataset(path_layermodel)

# %%
# Node property flow package - NPF
# =================================
#
# This package contains the information related to the aquifer properties used to calculate
# hydraulic conductance. This package replaces the Layer Property Flow (LPF),
# Block-Centered Flow (BCF), and Upstream Weighting (UPW) packages from previous MODFLOW versions.
#
# Hydraulic conductivity
# ----------------------
#
k = layermodel["k"]

npf = imod.mf6.NodePropertyFlow(
    icelltype=0,
    k=k,
    k33=k,
    variable_vertical_conductance=True,
    dewatered=True,
    perched=True,
    save_flows=True,
)

# %%
# Write
npf.to_netcdf(path_subsurface)
