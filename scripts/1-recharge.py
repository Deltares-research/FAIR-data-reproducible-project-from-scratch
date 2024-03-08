import imod
import numpy as np
import xarray as xr

# %%
path_meteorology = snakemake.input["path_meteorology"]
path_discretization = snakemake.input["path_discretization"]
path_recharge = snakemake.output["path_recharge"]

# %%
meteorology = xr.open_dataset(path_meteorology)
idomain = xr.open_dataset(path_discretization)["idomain"]

# %%
#
# Recharge
# ========
#
# This package is used to represent areally distributed recharge to the groundwater system.
# To calculate the recharge, the precipitation and evapotranspiration
# information from the KNMI website has been downloaded for the study area.
# This information is saved in netCDF files, which have been imported using the
# xarray function
# `xr.open_dataset <http://xarray.pydata.org/en/stable/generated/xarray.open_dataset.html#xarray.open_dataset>`_,
# slicing the area to the model's miminum and maximum dimensions.
#
# Note that the meteorological data has mm/d as unit and
# this has to be converted to m/d for Modflow 6.

xmin = 230_000.0
xmax = 257_000.0
ymin = 550_000.0
ymax = 567_000.0

pp = meteorology["precipitation"]
evt = meteorology["evapotranspiration"]

pp = pp.sel(x=slice(xmin, xmax), y=slice(ymax, ymin)) / 1000.0  # from mm/d to m/d
evt = evt.sel(x=slice(xmin, xmax), y=slice(ymax, ymin)) / 1000.0  # from mm/d to m/d

# %%
# Recharge - Steady state
# -----------------------
#
# For the steady state conditions of the model,
# the data from the period 2000 to 2009 was considered as reference.
# The initial information was sliced to this time period and averaged
# to obtain the a mean value grid. This process was done for both
# precipitation and evapotranspiration datasets.
#
# **Precipitation**
pp_ss = pp.sel(time=slice("2000-01-01", "2009-12-31"))
pp_ss_mean = pp_ss.mean(dim="time")

# %%
# **Evapotranspiration**
evt_ss = evt.sel(time=slice("2000-01-01", "2009-12-31"))
evt_ss_mean = evt_ss.mean(dim="time")


# %%
# For the recharge calculation, a first estimate
# is the difference between the precipitation and evapotranspiration values.

rch_ss = pp_ss_mean - evt_ss_mean


# %%
# Recharge - Transient
# --------------------
#
# The transient model will encompass the period from 2010 to 2015.
# The initial pp and evt datasets have been sliced to this time frame.

pp_trans = pp.sel(time=slice("2010-01-01", "2015-12-31"))
evt_trans = evt.sel(time=slice("2010-01-01", "2015-12-31"))

# %%
# As previously done, it is assumed that the recharge is equal
# to the difference between precipitation and evapotranspiration as a first estimate.
# Furthermore, the negative values found after doing this calculation have been
# replaced by zeros, as the recharge should not have a negative value.

rch_trans = pp_trans - evt_trans
rch_trans = rch_trans.where(rch_trans > 0, 0)  # check negative values

# %%
# The original information is on a daily step, so it is going to be
# resampled to a yearly step by using the xarray function
# `Dataset.resample <http://xarray.pydata.org/en/stable/generated/xarray.Dataset.resample.html#xarray.Dataset.resample>`_.

rch_trans_yr = rch_trans.resample(time="A", label="left").mean()
rch_trans_yr

# %%
# To create the final recharge for the transient simulation,
# the steady state information needs to be concatenated to the transient recharge data.
# The steady state simulation will be run for one second.
# This is achieved by using the numpy
# `Timedelta function <https://numpy.org/doc/stable/reference/arrays.datetime.html>`_,
# first creating a time delta of 1 second, which is assigned to the steady state recharge information.
# This dataset is then concatenated using the xarray function
# `xarray.concat <http://xarray.pydata.org/en/stable/generated/xarray.concat.html#xarray.concat>`_
# to the transient information and indicating that the dimension to join is "time".

starttime = "2009-12-31"

# Add first steady-state
timedelta = np.timedelta64(1, "s")  # 1 second duration for initial steady-state
starttime_steady = np.datetime64(starttime) - timedelta
rch_ss = rch_ss.assign_coords(time=starttime_steady)

rch_ss_trans = xr.concat([rch_ss, rch_trans_yr], dim="time")
rch_ss_trans

# %%
# The data obtained from KNMI has different grid dimensions
# than the one considered in this example. To fix this,
# imod-python includes the option
# :doc:`/api/generated/prepare/imod.prepare.Regridder`,
# which modifies the original grid dimensions to a different one.
# It is also possible to define the regridding method such as
# ``nearest``, ``multilinear``, ``mean``, among others.
# In this case, ``mean`` was selected and the 2d template (like_2d)
# was used as reference, as this is the geometry to be considered in the model.

like_2d = xr.full_like(idomain.isel(layer=0), 1)
like_2d

rch_ss_trans = imod.prepare.Regridder(method="mean").regrid(rch_ss_trans, like_2d)
rch_ss_trans

# %%
# The previously created recharge array is a 2D array
# that needs to be assigned to a 3D array. This is done using the xarray
# `DataArray.where <http://xarray.pydata.org/en/stable/generated/xarray.DataArray.where.html#xarray.DataArray.where>`_
# option, where the recharge values are applied to the cells where the
# idomain value is larger than zero (that is, the active cells) and for the uppermost
# active cell (indicated by the minimum layer number).

rch_total = rch_ss_trans.where(
    idomain["layer"] == idomain["layer"].where(idomain > 0).min("layer")
)
rch_total

# %%
# Finally, transposing the array dimensions using
# `DataArray.transpose <http://xarray.pydata.org/en/stable/generated/xarray.DataArray.transpose.html#xarray.DataArray.transpose>`_
# so they are in the correct order.

rch_total = rch_total.transpose("time", "layer", "y", "x")
rch_total

# %%
# Adding information to the RCH package
# --------------------------------------
#
# The information for the RCH package is added with the function
# :doc:`/api/generated/mf6/imod.mf6.Recharge`.
# It is required to insert the recharge flux rate, and it is optional
# to include the print_input, print_flows and save_flows information.

rch = imod.mf6.Recharge(rch_total)

# %%
rch.to_netcdf(path_recharge)