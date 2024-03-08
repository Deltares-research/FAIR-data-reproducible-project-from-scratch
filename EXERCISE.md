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
Systems. It features a ``README.md`` already describing the folder structure,
which you can further extend with a project description. Furthermore, it features
an ``AUTHORS.md`` file where all contributors to the project are credited. Also there is a ``LICENSE`` file describing the license.

> [!WARNING]
>
> The license added is very permissive, which we find works fine for most
> projects. However, you might want to change it to a stricter license. For
> example, for secret projects you probably want to add a propietary license.

All data that you start your project with is stored in the
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
is an essential part of software engineering these days and is also very useful
in our project work! The most common software for version control these days is
[git](https://git-scm.com/download/win). This has the advantage that it is very
well tested, documented, and a wealth of useful tools are available. For
example, some IDE's (e.g. VSCode) have a built-in git integration, allowing you
to version control your scripts from there. We'll run you through the basics and
will only skim the surface. Git allows a lot more things, for example easy
collaboration on a code base with colleagues.

# 3.1 Initializing git

We'll start off by initializing a git repository:

```powershell
git init
```

This will create a hidden ``.git`` folder, which contains the full history of
your scripts. 

# 3.2 Committing the initial state

You now have an empty version control system. Time to add some files to it!
First, let's see what files git can add to its version control system.
Type:

```powershell
git status
```

This will show you an overview of what files/folders git can add:

![git status initially](/docs/git_status_1.png)

To add all files, type:

```powershell
git add *
```

Check the status again:

```powershell
git status
```

This will show the files which are added:

![git status after add](/docs/git_status_2.png)

Now comes the most confusing part when learning git: adding files doesn't mean
they are safely stored yet in the version control system. For that we have to
commit:

```powershell
git commit -m "My initial commit"
```

After committing, files are added to the version control system. The reason why
git does this in two steps is that it allows you to orchestrate your commits
into logical steps for the history. In this case this is unnecessary, as we are
committing everything in one go, but is very useful in more complex situations
(you have to trust the millions of users git has these days on that for now...)

Type ``git status`` again and it will show you there's nothing to commit. You
have now succesfully safely stored your text files! Any change you make to the
files checked into git can now be tracked and reverted back to an old state.

# 3.3 Modifying a file and committing changes

Let's modify a file and commit the changes. Open the ``README.md`` in your
favorite editor and change the description. If you don't have inspiration what
to write, you can write "Reproducible workflow to run a groundwater model for
the Drentse Hondsrug". Save the file, and type ``git status`` again to confirm
git noticed changes to the file.

![git status after changes](/docs/git_status_3.png)

Note that git also tips you with some commands you could use next: You can
either ``add`` the README to store changes, or ``restore`` it to its last
committed state.

Next, let's review the changes you made to the file exactly. Type:

```powershell
git diff
```

This will show you the changes you made to the text:

![git diff](/docs/git_diff.png)

If you're satisfied with these changes, we can add them:

```powershell
git add README.md
```

and commit them: 

```powershell
git commit -m "Modify project description"
```

# 3.4 Checking version history

We can keep track of our version history by typing:

```powershell
git log
```

This will print you the two commits you made with their commit messages. Note
that it therefore is important to write short but descriptive commit messages,
so you can more easily retrace your steps! Commit messages like "Update" are too
generic to be of any use.

# 3.5 Excluding files from source version control

Note that your folder also contains a ``.gitignore`` file. This was included in
the cookiecutter project template. It can be opened and edited in any text
editor. Open it in your favorite text editor, and you will see certain folders
and file extensions being listed here. These files will be ignored by git, and
thus not kept in version control. This is useful: We do not want to store all
our files, as bulky files can easily clog the version control system, and don't
have to be stored. We already added a bunch of common files you are very likely
not to want to add to your version control system in here. For example, you do
not want your ``.pixi`` folder, containing 2GB of python installation specific
to your machine, checked in git: the ``pixi.toml`` and ``pixi.lock`` file are
enough to recreate the ``.pixi`` folder. Therefore the ``.pixi`` folder is
included in the ``.gitignore`` file. In general, it is best to not commit large
files to git. Regular git also is not very useful to work with binary data. In
that case, you are better off using git-lfs or DVC. We'll explain how to use DVC
in a later stage of this exercise.

# 3.6 Adding scripts to repository

Finally, we'll add the Python scripts, which are our workflow steps. These are
already prepared in the folder ``scripts``. The scripts are named with a prefix
number indicating in which folder under ``src`` they are supposed to be put. For
example, ``0-download-data.py`` should be moved to the folder ``src/0-setup``.
Copy all scripts to their respective folder.

If everything went well ``git status`` will list the following files:

![git status scripts](/docs/git_status_scripts.png)

Add and commit these files the same way you did this before:

```powershell
git add *
git commit -m "Added scripts to repository"
```

# 4 Setting up the workflow: Snakemake 

We have a collection of scripts, which are depending on each other. For example
``0-download-data.py`` should be called before calling ``1-surface-water.py``,
as downloaded data is required to schematize the surface water system. Snakemake
takes care of this. How it works is by checking for each step in the workflow
what data comes in and what data comes out. This has to be specified explicitly
by the user. For example, we have to tell snakemake that the file ``river.nc``
is output of the script ``0-download-data.py`` and input to the script
``1-surface-water.py``. Snakemake then will deduce by itsself that
``0-download-data.py`` has to be called before ``1-surface-water.py``. This
might seem underwhelming for such a trivial situation, but it gets very useful
in more complex situations, as snakemake deduces the dependence of steps and
order of computation by itsself.

# 4.1 Create a snakefile

To start configuring our snakemake workflow, start off by creating a file named
``snakefile``. By default, snakemake will look for a file named ``snakefile`` as
its configuration file. Snakemake defines its individual steps as "rules". Let's
add our first rule to the ``snakefile``. Open up your favorite editor, and copy
the following rule in your ``snakefile``:

```
rule download_data:
    output:
        path_layermodel = "data/1-external/layermodel.nc",
        path_starting_heads = "data/1-external/starting_heads.nc",
        path_meteorology = "data/1-external/meteorology.nc",
        path_drainage = "data/1-external/drainage.nc",
        path_river = "data/1-external/river.nc",
    script:
        "src/0-setup/0-download-data.py"
```

This rule will call the script ``src/0-setup/0-download-data.py`` and checks if
it produced the files: ``layermodel.nc``, ``starting_heads.nc``,
``meteorology.nc``, ``drainage.nc``, ``river.nc``.

Now run:

```powershell
snakemake -c1
```

This will run the snakemake workflow. The option ``-c1`` is shorthand for
``--cores 1`` and will thus tell snakemake to use only one core of your machine.
This is enough, as we have not defined any independent steps which can be run in
parallel.

# 4.2 Adding a second step

Let's define the second step to our workflow. We'll call the script
``src/1-prepare/1-discretization.py`` which creates the model's spatial
discretization for Modflow 6, based on hydrogeological layers provided in
``layermodel.nc``. We'll add the rule ``discretization`` above the
``download_data`` rule, as snakemake by default will look at the first rule to
run it (and all its dependencies.)

```
rule discretization:
    input:
        path_layermodel = "data/1-external/layermodel.nc",
    output:
        path_discretization = "data/2-interim/discretization.nc",
    script:
        "src/1-prepare/1-discretization.py"

rule download_data:
    output:
        path_layermodel = "data/1-external/layermodel.nc",
        path_starting_heads = "data/1-external/starting_heads.nc",
        path_meteorology = "data/1-external/meteorology.nc",
        path_drainage = "data/1-external/drainage.nc",
        path_river = "data/1-external/river.nc",
    script:
        "src/0-setup/0-download-data.py"
```

Run again in powershell:

```powershell
snakemake -c1
```

# 4.2 Finish your snakefile

<details>
  <summary> If you did everything correct, your Snakefile will look as follows: (<i>click to expand</i>)</summary>
  <!-- have to be followed by an empty line! -->

```
rule plot_heads:
    input:
        path_head_nc = "data/5-visualization/groundwater_heads.nc"
    output:
        path_figure = "reports/figures/groundwater_heads.png"
    script:
        "src/5-visualize/5-plot.py"

rule post_process:
    input:
        path_hds = "data/4-output/GWF.hds",
        path_grb = "data/4-output/dis.dis.grb",
    output:
        path_head_nc = "data/5-visualization/groundwater_heads.nc"
    script:
        "src/4-analyze/4-post-process.py"

rule run_model:
    input:
        path_model = "data/3-input/mfsim.nam"
    output:
        path_hds = "data/4-output/GWF.hds",
        path_grb = "data/4-output/dis.dis.grb",
    shell:
        "cd data\\3-input && call ..\\..\\bin\\mf6.exe . && move GWF\\GWF.hds ..\\4-output\\GWF.hds && move GWF\\dis.dis.grb ..\\4-output\\dis.dis.grb"

rule build_model:
    input:
        path_discretization = "data/2-interim/discretization.nc",
        path_drn_pkg = "data/2-interim/drn_pkg.nc",
        path_riv_pkg = "data/2-interim/riv_pkg.nc",
        path_recharge =  "data/2-interim/recharge.nc",
        path_ic = "data/2-interim/ic.nc",
        path_chd = "data/2-interim/chd.nc",
        path_subsurface = "data/2-interim/subsurface.nc",
    output:
        path_model = "data/3-input/mfsim.nam"
    script:
        "src/2-build/2-build-model.py"

rule surface_water:
    input:
        path_drainage = "data/1-external/drainage.nc",
        path_river = "data/1-external/river.nc",
    output:
        path_drn_pkg = "data/2-interim/drn_pkg.nc",
        path_riv_pkg = "data/2-interim/riv_pkg.nc",
    script:
        "src/1-prepare/1-surface-water.py"

rule recharge:
    input:
        path_meteorology = "data/1-external/meteorology.nc",
        path_discretization = "data/2-interim/discretization.nc",
    output:
        path_recharge = "data/2-interim/recharge.nc",
    script:
        "src/1-prepare/1-recharge.py"

rule initial_condition:
    input:
        path_starting_heads = "data/1-external/starting_heads.nc",
        path_discretization = "data/2-interim/discretization.nc",
    output:
        path_ic = "data/2-interim/ic.nc",
        path_chd = "data/2-interim/chd.nc",
    script:
        "src/1-prepare/1-initial-condition.py"

rule subsurface:
    input:
        path_layermodel = "data/1-external/layermodel.nc",
    output:
        path_subsurface = "data/2-interim/subsurface.nc",
    script:
        "src/1-prepare/1-subsurface.py"

rule discretization:
    input:
        path_layermodel = "data/1-external/layermodel.nc",
    output:
        path_discretization = "data/2-interim/discretization.nc",
    script:
        "src/1-prepare/1-discretization.py"

rule download_data:
    output:
        path_layermodel = "data/1-external/layermodel.nc",
        path_starting_heads = "data/1-external/starting_heads.nc",
        path_meteorology = "data/1-external/meteorology.nc",
        path_drainage = "data/1-external/drainage.nc",
        path_river = "data/1-external/river.nc",
    script:
        "src/0-setup/0-download-data.py"

```


