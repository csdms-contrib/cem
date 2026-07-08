[![Build Status](https://travis-ci.org/csdms-contrib/cem.svg?branch=master)](https://travis-ci.org/csdms-contrib/cem)

The Coastline Evolution Model (CEM)
===================================

The Coastline Evolution Model (or CEM) is a one-contour line model that focuses
on sandy, wave-dominated shoreline evolution, simulates the plan-view evolution
of a coastline due to gradients in alongshore sediment transport. A unique
aspect of CEM, by dividing the plan-view domain into a 2-dimensional cell array,
is its ability to process an arbitrarily sinuous shoreline, allowing the
simulation of complex shoreline features including spits and capes.  The model
is exploratory in nature, designed to simulate large-scale (10^3 to 10^6 m) and
long-term (10^2 to 10^5 yr) shoreline evolution.

Build and install CEM
---------------------

CEM can be built from source on Linux, macOS, and Windows.
Instructions are given below.

**Prerequisites:**
* A C compiler
* CMake

A convenient way to build and install CEM is within a
[virtual environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
set up through the [conda](https://docs.conda.io/projects/conda/en/latest/index.html) package manager.
After cloning or downloading this repository,
change into the repository directory
and set up a conda environment with the included environment file:
```sh
conda env create --file environment.yml
```

To build CEM from source with CMake, run
```sh
cmake -B _build --fresh
cmake --build _build
```
The compiled model files live in the directory `_build`.

Run the tests with
```sh
ctest -V --test-dir _build
```

Install the model with
```sh
cmake --install _build --prefix PATH-TO-INSTALLATION
```
where `PATH-TO-INSTALLATION` is the base directory
in which to install the model (`/usr/local` is the default).
When using a conda environment,
use the `$CONDA_PREFIX` environment variable.

Build the old version of CEM
----------------------------

Build the old version of CEM,

    > gcc main.c -o cem.exe

To run this version of CEM, it's probably best to run it in another directory
as it generates lots of output files,

    > mkdir _run
    > cd _run
    > cem.exe
