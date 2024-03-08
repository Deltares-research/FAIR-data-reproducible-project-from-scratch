import matplotlib.pyplot as plt
import xarray as xr

path_head_nc = snakemake.input["path_head_nc"]
path_figure = snakemake.output["path_figure"]

# %%
heads = xr.open_dataarray(path_head_nc)

# %% Plot
fig, ax = plt.subplots()
heads.sel(layer=slice(3, 5)).mean(dim="layer").isel(time=3).plot(ax=ax)

# %% Write
plt.savefig(path_figure)
