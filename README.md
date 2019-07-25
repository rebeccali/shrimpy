# SHRIMP
Repository for SHRIMP robot dynamics simulation.

Maintained by Rebecca Li (robot at seas.upenn.edu) and Allen Yu (fangzhyu at seas.upenn.edu)


# Repository Structure

The new simulation code is in `src`. You can quickly run a simulation by running `runSim.py`. It is connected to many of the other files in the `src` folder. For this simulation, you will need the conda environment which you can install as described in the `Python Installation`section.


The second simulator is `ShrimpSimulator`, which is simply a python notebook and self explanatory. For `ShrimpSimulator`, you will need the conda environment which you can install as described in the `Python Installation` section.


The old simulation code is contained within `piccolissimo`. Each simulation is slighlty different:

* `UNO` is the simulator in active development for the SHIRMP project.
* `generic_sim` is for a generic simulator, but is difficult to use.
* `simple_sim` is a simple simulator, but is difficult to use.

## Running a native python simulation:
You can either run or test the code. To run the code, simply type

`python src/runSim.py`

During development, be sure to also run

`python src/test.py`

This will be eventually integrated into CI

## Running `ShrimpSimulator` Notebook :

0. Make sure Jupyter notebook is installed.  This comes with Anaconda, but otherwise follow installation instructions from https://jupyter.org/.
1. Run `jupyter notebook` in this directory.
2. From the browser interface, open up `src/ShrimpSimulator.ipynb`.
3. Run all cells.

## Running an UNO simulation:

Currently, there are two ways to run simulations

* `RunBasicSimulator.m` runs a basic half-second simulation and generates a bunch of plots to file. It is very simple, and should only be used for regression testing.
* `RunOptimizeBetaD.m` runs a series of simulations in order to calculate the shape of the body stabilizers. It is not actively maintained.


## Python Installation
The python aspects requires python.

1. Install Miniconda or Anaconda https://docs.conda.io/en/latest/miniconda.html
2. Install the python environment `conda env create -f environment.yml -n`.
3. Test that the environment exists by opening a new terminal and running `source activate shrimp`. This is the new python environment with libraries preloaded.
4. When running any of the `ShrimpSimulator` files, use this environment.
5. During development, be sure to update `environment.yml` if you install any new packages.

## The `src/matlab` folder
An old attempt at matlab simulation.
TODO(Allen): documentation

## Tasks to Do:
* Set up CI
* Test out
