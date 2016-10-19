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

Build CEM
---------

First create a build directory and then use [cmake](http://cmake.org)
to generate a makefile,

    > mkdir _build
    > cd _build
    > cmake ..

Use the makefile to compile CEM,

    > make

Continuous Integration
----------------------

Continuous integration is handled by
[Travis-CI](https://travis-ci.org/csdms-contrib/cem). Every push to this
GitHub repository triggers a new build on Travis. Travis builds CEM on
both Mac and Linux and runs a simple smoke test to see that things are
working. The `.travis.yaml` file controls how Travis builds and runs
tests.

Build the old version of CEM
----------------------------

Build the old version of CEM,

    > gcc main.c -o cem.exe

To run this version of CEM, it's probably best to run it in another directory
as it generates lots of output files,

    > mkdir _run
    > cd _run
    > cem.exe

