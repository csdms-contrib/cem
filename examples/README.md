CEM Examples
============

This folder contains example files for the coupled Waves and CEM
models. The `cem` program requires two input files:

*  `waves.txt` defines parameters for the waves component
*  `cem.txt` defines parameters for the CEM component

The provided files in this folder contain some default parameters to
get you going. After installing both components, you can run the
example with the following command:

    $ cem waves.txt cem.txt output.csv

This will create an `output.csv` file that will contain snapshots of
the evolving delta at regular intervals.

If you're wondering what the parameters in the input files are, have
a look at the template input files (`.tmpl`) and parameter definition
files (`parameters.yaml`) in `.bmi-waves` and `.bmi-cem`.
