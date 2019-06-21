# SHRIMP
Repository for SHRIMP robot dynamics simulation.

Maintained by Allen Yu (fangzhyu at seas.upenn.edu) and Rebecca Li (robot at seas.upenn.edu).


# Repository Structure

The simulation code is contained within `piccolissimo`. Each simulation is slighlty different:

* `UNO` is the simulator in active development for the SHIRMP project.
* `generic_sim` is for a generic simulator, but is difficult to use.
* `simple_sim` is a simple simulator, but is difficult to use.

## Running an UNO simulation:

Currently, there are two ways to run simulations

* `RunBasicSimulator.m` runs a basic half-second simulation and generates a bunch of plots to file. It is very simple, and should only be used for regression testing.
* `RunOptimizeBetaD.m` runs a series of simulations in order to calculate the shape of the body stabilizers. It is not actively maintained.


## Tasks to Do:
* Remove globals
* Make a better simulation visualization
* Clean up configuration files
* Write documentation on using configuration files
