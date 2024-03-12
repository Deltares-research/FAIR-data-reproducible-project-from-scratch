# Introduction

This exercise is the second part of the FAIR data reproducible projects pizza
course, where we create a reproducible project from scratch. It is the follow-up
to [the first exercise, which can be found
here,](https://github.com/Deltares-research/FAIR-data-example-project) where we
run a reproducible project on our local machine.

# Exercise

In the exercise we'll create a reproducible from scratch. The starting point is
a set of scripts, which already have been prepared for you and can be found in
the [scripts folder](scripts), and data which will be downloaded from the
internet. Furthermore, Modflow 6 the executable has to be downloaded manually,
which will be explained in the exercise.

[Link to the exercise materials.](EXERCISE.md)

# Context

The context is not very important for this exercise, but for those interested
some background will be described here. The workflow will run a groundwater flow
model for the Drentse Hondsrug.

The "Hondsrug" a Dutch ridge of sand over a range of 70 km. It is the only
*geopark* in The Netherlands and part of Natura 2000, an European network of
protected areas.

Scripts were created to prepare and analyze a groundwater flow model. The model
is based on a test model in the iMOD5 test bench, which is loosely based on an
old version of a MIPWA model. This is the regional groundwater of the North-East
of the Netherlands, used and maintained by several governmental agencies. [More
info about this regional model here.](https://nhi.nu/modellen/mipwa/) The
applied model code is Modflow 6, which is developed by the USGS. [Find more info
about this model code here.](https://github.com/MODFLOW-USGS/modflow6)

![foto_hondsrug](docs/hondsrug.jpg)