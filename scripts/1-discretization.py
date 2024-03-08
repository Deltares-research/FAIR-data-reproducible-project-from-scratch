# %%
import imod
import xarray as xr

# %%
# Path management
path_layermodel = snakemake.input["path_layermodel"]
path_discretization = snakemake.output["path_discretization"]

# %%
# This package allows specifying a regular MODFLOW grid. This grid is assumed
# to be rectangular horizontally, but can be distorted vertically.
#
# Load data
# ---------
#
# We'll load the data from the examples that come with this package.
layermodel = xr.open_dataset(path_layermodel)

# Make sure that the idomain is provided as integers
idomain = layermodel["idomain"].astype(int)

# We only need to provide the data for the top as a 2D array. Modflow 6 will
# compare the top against the uppermost active bottom cell.
top = layermodel["top"].max(dim="layer")

bot = layermodel["bottom"]

# %% Create discretization object
dis = imod.mf6.StructuredDiscretization(
    top=top, bottom=bot, idomain=idomain
)

dis.to_netcdf(path_discretization)
