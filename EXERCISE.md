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

## 2.1 Prepare python environment: Pixi

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
Inspect your folder's contents in TotalCommander/Windows Explorer.
Alternatively, to inspect folder contents, you can print files and folders in
your shell session by calling in Powershell:

```powershell
dir
```

You can see the ``pixi init`` command created some text files, the important one
now being the ``pixi.toml`` file, which is a configuration file, containing all
the important settings to create a pixi python environment. We'll see how this
file gets extended later.

Next, add cookiecutter to your project dependencies:

```powershell
pixi add cookiecutter jinja2-time
```

> [!NOTE]
>
> At the time of writing (2024-03-06), cookiecutter didn't automatically install
> a dependency ``jinja2-time``. Therefore this package has to be added manually.

Inspect your folder again. Pixi created a hidden folder ``.pixi`` and a
``pixi.lock`` file, containing the python environment and text representation of
the exact state of the python environment contents. 

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

Now activate your pixi environment in a shell session:

```powershell
pixi shell
```

## 2.2 Apply project template: Cookiecutter

### 2.2.1 Running cookiecutter

Next we'll apply the project template created by the Groundwater Management
Department, run:

```powershell
cookiecutter gl:deltares/imod/cookiecutter-reproducible-project
```

This will ask you to fill in some details about the project and consequently
creates a folder structure.

### 2.2.2 Moving pixi environment into project safely

We'll inspect the folder in detail soon, but first we have to deal with one
minor inconvenience, namely that it is the most convenient to have the pixi
files in the project root folder (i.e. one folder down). Therefore, move the
``pixi.toml`` and ``pixi.lock`` file to the project folder. You can do this
manually in your explorer, or do it from powershell:

```powershell
mv pixi.* <your_project_name>/
```

> [!NOTE]
>
> The asterisk (*) acts as a wildcard, so all files starting with "pixi." are
> matched (and moved), in this case ``pixi.lock`` and ``pixi.toml``.

Let's declutter some more by removing the pixi environment, we'll recreate it
later!

First exit your pixi session:

```powershell
exit
```

Then remove the ``.pixi`` folder:

```powershell
rm .pixi -Force -Recurse
```

> [!TIP]
>
> A python environment consists of a lot of files, easily over 40K, and easily
> becomes quite large, because of some common dependencies: most notorious is
> the ``mkl`` dependency. Therefore make sure to always permanently delete your
> ``.pixi`` environment, instead of it being moved to the Recylce Bin. In the
> Windows Explorer/Total Commander, this is done with shortcut key combination
> SHIFT+DELETE.

Move into your project folder:

```powershell
cd <your_project_name>
```

Before we continue to the next step, let's recreate our python environment
again:

```powershell
pixi install
```

Notice that it takes very little time for pixi to create your python environment
again! Now activate your pixi environment in a shell session:

```powershell
pixi shell
```

### 2.2.3 Project template folder structure

Inspect the folder again by calling:

```powershell
dir
```

You'll see the following structure (also conveniently described in the
``README.md``):

    .
    ├── AUTHORS.md
    ├── LICENSE
    ├── README.md
    ├── bin                 <- Your compiled model code can be stored here (not tracked by git)
    ├── config              <- Configuration files, e.g., for doxygen or for your model if needed
    ├── data                
    │   ├── 1-external      <- Data external to the project.
    │   ├── 2-interim       <- Intermediate data that has been altered.
    │   ├── 3-input         <- The processed data sets, ready for modeling.
    │   ├── 4-output        <- Data dump from the model.
    │   └── 5-visualization <- Post-processed data, ready for visualisation.
    ├── docs                <- Documentation, e.g., doxygen or scientific papers (not tracked by git)
    ├── notebooks           <- Jupyter notebooks
    ├── reports             <- For a manuscript source, e.g., LaTeX, Markdown, etc., or any project reports
    │   └── figures         <- Figures for the manuscript or reports
    └── src                 <- Source code for this project
        ├── 0-setup         <- Install necessary software, dependencies, pull other git projects, etc.
        ├── 1-prepare       <- Scripts and programs to process data, from 1-external to 2-interim.
        ├── 2-build         <- Scripts to create model specific inputm from 2-interim to 3-input. 
        ├── 3-model         <- Scripts to run model and convert or compress model results, from 3-input to 4-output.
        ├── 4-analyze       <- Scripts to post-process model results, from 4-output to 5-visualization.
        └── 5-visualize     <- Scripts for visualisation of your results, from 5-visualization to ./report/figures.

This project template is commonly applied in the unit Subsurface and Groundwater
Systems. All data that you start your project with is stored in the
``data/1-external`` folder, these are usually the files you received from a
client, or downloaded somewhere. In case of model update, you could treat the
model before the update as "external data". This data is the starting point of
your workflow. Usually external data has to be reworked and cleaned up in order
to lead to meaningful results, this leads to pre-processed data is stored in
``data/2-interim``. Most model codes have their own specific files (often not
interoperable), these are stored ``data/3-input``. Model output data is stored
in ``data/4-output``. Finally, model output has to be post-processed for
plotting, this for example can be converting model-specific formats to more
interoperable file formats, converting a 3D grid to VTK blocks for 3D plotting,
or aggregating the data into a timeseries for a line plot. This data is stored
in ``data/5-visualization``. Figures created from this post-processed data are
stored in ``reports/figures``.

> [!TIP]
>
> Some model codes can only write model output in the input folder, in that case
> it is wise to add a script/command to move output files from ``data/3-input``
> to ``data/4-output``

# 3 Source version control: git

Version control means keeping track of all the changes made to the project. It
is an essential part of software engineering these days, but is also very useful
in our project work! The most common software for version control these days is
git. This has the advantage that it is very well tested and documented.

# 3.1 Initializing git

We'll start off by initializing a git repository:

```powershell
git init
```

This will create a hidden ``.git`` folder, which contains the full history of
your scripts. Note that your folder also contains a ``.gitignore`` file. This
was created in the call to cookiecutter template. 

This can be opened and edited in any text editor. 

# 3.2 Adding scripts

Up next, we'll add some Python scripts, which are our workflow steps. If you 

```powershell

```
