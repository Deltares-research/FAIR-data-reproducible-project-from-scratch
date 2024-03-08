import imod
import numpy as np
import scipy
import xarray as xr

# %%
# Path management
path_discretization = snakemake.input["path_discretization"]
path_starting_heads = snakemake.input["path_starting_heads"]
path_initial = snakemake.output["path_ic"]
path_constant_head = snakemake.output["path_chd"]

# %%
starting_heads = xr.open_dataset(path_starting_heads)
idomain = xr.open_dataset(path_discretization)["idomain"]

# %%
# Initial conditions package - IC
# ================================
#
# This package reads the starting heads for a simulation.
#
# Starting heads interpolation
# ----------------------------
#
# The starting heads to be used in this model are based on the interpolation of
# x-y head measurements, which were interpolated on a larger area.  This
# example was created in this example --insert-link-here--
#
# The heads were interpolated on a larger area, therefore these have to be
# clipped first

interpolated_head_larger = starting_heads["head"]

xmin = 237_500.0
xmax = 250_000.0
ymin = 559_000.0
ymax = 564_000.0

interpolated_head = interpolated_head_larger.sel(
    x=slice(xmin, xmax), y=slice(ymax, ymin)
)

# %%
# The final step is to assign the 2D heads interpolation to all the
# model layers (as a reference value) by using the xarray tool
# `xarray.full_like <http://xarray.pydata.org/en/stable/generated/xarray.full_like.html#xarray.full_like>`_.
# The 3d idomain array is used as reference for the geometry and then
# its original values are replaced by NaNs.
# This array is combined with the interpolated_head array using the xarray
# `DataArray.combine_first <http://xarray.pydata.org/en/stable/generated/xarray.DataArray.combine_first.html#xarray.DataArray.combine_first>`_
# option.
# The final result is an starting_heads xarray where all layers have the 2d interpolated_head information.

# Assign interpolated head values to all the model layers
like_3d = xr.full_like(idomain, np.nan, dtype=float)
starting_head = like_3d.combine_first(interpolated_head)
# Consequently ensure no data is specified in inactive cells:
starting_head = starting_head.where(idomain == 1)

starting_head

# %%
# Adding information to the IC package
# ------------------------------------
#
# The function for indicating the initial conditions is
# :doc:`/api/generated/mf6/imod.mf6.InitialConditions`.
# It is necessary to indicate the value(s) to be considered as the initial
# (starting) head of the simulation.
# In this case, this value is equal to the previously created starting_head array.

ic = imod.mf6.InitialConditions(start=starting_head)

# %%

# %%
# Constant head package - CHD
# ===========================
#
# This package allows to indicate if the head varies with time,
# if it is constant or if it is inactive.
#
# Constant head edge
# -------------------
#
# The previously interpolated starting_head array will be used to define
# the constant head value which will be used along the model boundaries.
# A function is defined to indicate the location of the outer edge
# (returning a boolean array).


def outer_edge(da):
    data = da.copy()
    from_edge = scipy.ndimage.binary_erosion(data)
    is_edge = (data == 1) & (from_edge == 0)
    return is_edge.astype(bool)


# %%
# For the next calculations, it is necessary to create a template array
# which can be used for assigning the corresponding geometry to other arrays.
# In this case, a 2d template is created based on the idomain layer information
# and filled with ones.

like_2d = xr.full_like(idomain.isel(layer=0), 1)
like_2d

# %%
# Using the previously created function and the 2d template,
# the outer edge is defined for this example.

edge = outer_edge(xr.full_like(like_2d.drop_vars("layer"), 1))

# %%
# Adding information to the CHD package
# --------------------------------------
#
# To add the information to the CHD package within the gwf_model variable, the
# :doc:`/api/generated/mf6/imod.mf6.ConstantHead`.
# function is used.
# The required information is the head array for this boundary condition.
# In this example, the starting_head array is selected where the idomain is > 0 (active)
# and it is located in the edge array.
#
# It is also possible (and optional) to indicate if the CHD information will be written
# to the listing file after it is read (print_input), if the constant head flow rates will
# be printed to the listing file for every stress period time step
# in which “BUDGET PRINT” is specified in Output Control (print_flows)
# and if the constant head flow terms will be written to the file
# specified with “BUDGET FILEOUT” in Output Control (save_flows).
# By default, these three options are set to False.

chd = imod.mf6.ConstantHead(
    starting_head.where((idomain > 0) & edge),
    print_input=False,
    print_flows=True,
    save_flows=True,
)

# %%
# Write
ic.to_netcdf(path_initial)
chd.to_netcdf(path_constant_head)

