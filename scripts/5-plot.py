import matplotlib.pyplot as plt

plt.switch_backend("Agg")
import contextily as ctx
import imod
import numpy as np
import xarray as xr

path_head_nc = snakemake.input["path_head_nc"]
path_figure = snakemake.output["path_figure"]

# %% Open data
heads = xr.open_dataarray(path_head_nc)

# %% Select data to plot
data_to_plot = heads.sel(layer=slice(3, 5)).mean(dim="layer").isel(time=-1)

# %% Plot
fig, ax = plt.subplots()

# Access background maps
background_map = ctx.providers["OpenStreetMap"]["Mapnik"]

colors = "viridis"
levels = np.arange(-2.5, 18.5, 2.5)

imod.visualize.plot_map(data_to_plot, colors, levels, basemap=background_map, fig=fig, ax=ax)

ax.set_title("Groundwater heads (m +NAP)")

# %% Write
plt.savefig(path_figure)
