# SHRIMP
Repository for SHRIMP robot dynamics simulation.

Maintained by Allen Yu (fangzhyu at seas.upenn.edu) and Rebecca Li (robot at seas.upenn.edu).


# Repository Structure

The new simulation code is in `src`. The latest simulator is `RageSim`, which is simply a python notebook and self explanatory. For `RageSim`, you will need the conda environment which you can install as described in the `Python Installation section`


The old simulation code is contained within `piccolissimo`. Each simulation is slighlty different:

* `UNO` is the simulator in active development for the SHIRMP project.
* `generic_sim` is for a generic simulator, but is difficult to use.
* `simple_sim` is a simple simulator, but is difficult to use.

## Running an UNO simulation:

Currently, there are two ways to run simulations

* `RunBasicSimulator.m` runs a basic half-second simulation and generates a bunch of plots to file. It is very simple, and should only be used for regression testing.
* `RunOptimizeBetaD.m` runs a series of simulations in order to calculate the shape of the body stabilizers. It is not actively maintained.


## Installing Python for `RageSim`
RageSim requires python.

1. Install Miniconda or Anaconda https://docs.conda.io/en/latest/miniconda.html
2. Install the python environment `conda env create -f environment.yml -n`.
3. Test that the environment exists by opening a new terminal and running `source activate shrimp`. This is the new python environment with libraries preloaded.
4. When running any of the `RageSim` files, use this environment.
5. During development, be sure to update `environment.yml` if you install any new packages.

## Tasks to Do:
* Remove globals
* Make a better simulation visualization
* Clean up configuration files
* Write documentation on using configuration files

