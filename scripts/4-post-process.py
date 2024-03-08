import imod

# Path management
path_hds = snakemake.input["path_hds"]
path_grb = snakemake.input["path_grb"]
path_head_nc = snakemake.output["path_head_nc"]

# Open heads
heads = imod.mf6.out.open_hds(path_hds, path_grb)

# Write heads
heads.to_netcdf(path_head_nc)