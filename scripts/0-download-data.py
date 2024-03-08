import imod

# Path management
path_layermodel = snakemake.output["path_layermodel"]
path_initial = snakemake.output["path_starting_heads"]
path_meteorology = snakemake.output["path_meteorology"]
path_drainage = snakemake.output["path_drainage"]
path_river = snakemake.output["path_river"]

# Download imod example data with Pooch
layermodel = imod.data.hondsrug_layermodel()
initial = imod.data.hondsrug_initial()
meteorology = imod.data.hondsrug_meteorology()
drainage = imod.data.hondsrug_drainage()
river = imod.data.hondsrug_river()

# Save files
layermodel.to_netcdf(path_layermodel)
initial.to_netcdf(path_initial)
meteorology.to_netcdf(path_meteorology)
drainage.to_netcdf(path_drainage)
river.to_netcdf(path_river)
