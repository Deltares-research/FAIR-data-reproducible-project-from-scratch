# Exercise 2: Creating a reproducible from scratch

This excercise is the second card of the reproducible modelling pizza course.
The first part of the course demonstrated how to work with an existing
reproducible model and [can be found
here](https://github.com/Deltares-research/FAIR-data-example-project).

In this second part of the course, we'll learn how to create a reproducible
workflow from scratch. 

## Requirements

We expect you to have installed:

- [pixi](https://pixi.sh/latest/)
- [git](https://git-scm.com/download/win)

# 1. Create directory

First create an empty folder somewhere on your machine. In general we advise to
do this outside the OneDrive folder, as that has the following downsides:
    
- The folder contains spaces (at least on Deltares laptops), this easily leads
  to mistakes as this often requires extra user input (e.g. put the path between
  quotes). Furthermore often software contains bugs or doesn't even support
  spaces in paths.
- A python installation will be created, which creates lots of files (>40k
  easily). This can easily clog your OneDrive synchroniztion
- It is not necessary, our project will be reproducible!

Open a powershell/cmd session and type:

```powershell
mkdir c:/users/<your_name>/<path>/<to>/<folder>
```

Next, navigate into the folder:

```powershell
cd c:/users/<your_name>/<path>/<to>/<folder>
```

# 2. Initialize project

## 2.1 Prepare python environment

This workflow requires the following dependencies:

- [iMOD Python](https://deltares.github.io/imod-python/): Generate and
  process Modflow 6 models
- [DVC](https://dvc.org/): Data Version Control
- [cookiecutter](https://cookiecutter.readthedocs.io/en/stable/): Apply project templates
- [snakemake](https://snakemake.readthedocs.io/en/stable/): Workflow manager

First initialize pixi **in your project folder**. We require the use of two
channels. Most packages are located on ``conda-forge``, but Snakemake is only
published on ``bioconda``, therefore we specify this channel as well.

```powershell
pixi init --channel conda-forge --channel bioconda
```

Next, add cookiecutter to your project dependencies:

```powershell
pixi add cookiecutter
```

Inspect your folder, pixi created a hidden folder ``.pixi``, containing the
python environment.

We'll add DVC. During the creation of this course, we found that dvc tended to
install an older version by default, therefore it is best to force installing
a later version:

```powershell
pixi add "dvc>=3.48.2"
```

Next, add iMOD Python. To save our iMOD developer colleagues some work for any
future breaking changes, we'll force you to install the latest version at time
of writing this material (2023-03-05).

```powershell
pixi add "imod=0.15.3"
```

Finally, add Snakemake. This will be installed from the bioconda channel.

```powershell
pixi add snakemake
```

Now let's lose some work! Let's "accidentily" remove our pixi environment:

```powershell
rm .pixi -Force -Recurse
```

Oh no! We've lost our python installation! Normally this results in long waiting
times to install again. Luckily pixi makes reinstalling very fast. Let's install
our pixi environment again:

```powershell
pixi install
```

Now activate your pixi environment in a shell session:

```powershell
pixi shell
```

## 2.2 Apply project template

Next we'll apply the project template created by the Groundwater Management
Department, run:

```powershell
cookiecutter gl:deltares/imod/cookiecutter-reproducible-project
```

```powershell
git init
```

```powershell

```
