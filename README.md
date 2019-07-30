# SHRIMP
Repository for SHRIMP robot dynamics simulation.

Maintained by Rebecca Li (robot at seas.upenn.edu) and Allen Yu (fangzhyu at seas.upenn.edu)

# QuickStart

## Python Installation
The python aspects requires python.

1. Install Miniconda or Anaconda https://docs.conda.io/en/latest/miniconda.html
2. Install the python environment `conda env create -f environment.yml -n`.
3. Test that the environment exists by opening a new terminal and running `source activate shrimp`. This is the new python environment with libraries preloaded.
4. When running any of the `ShrimpSimulator` files, use this environment.
5. During development, be sure to update `environment.yml` if you install any new packages.


Note that should you ever need to update the conda environment, after updating run the following command:

`conda env export > environment,yml`

## Run a test:

`python src/runSim.py`

Or to go through the whole codebase:

`python src/test.py`


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

This will be eventually integrated into CI. This code is heavily a port from the original `UNO` simulation, but has different notation, which can be confusing.

## Running `ShrimpSimulator` Notebook :

0. Make sure Jupyter notebook is installed.  This comes with Anaconda, but otherwise follow installation instructions from https://jupyter.org/.
1. Run `jupyter notebook` in this directory.
2. From the browser interface, open up `src/ShrimpSimulator.ipynb`.
3. Run all cells.

## Running an UNO simulation:

Currently, there are two ways to run simulations

* `RunBasicSimulator.m` runs a basic half-second simulation and generates a bunch of plots to file. It is very simple, and should only be used for regression testing.
* `RunOptimizeBetaD.m` runs a series of simulations in order to calculate the shape of the body stabilizers. It is not actively maintained.

### Tips on reading through the original UNO simulation
TODO: this, as well as discussing the difference in notation, such as `Rr_f` vs `rot_p2f`


## The `src/matlab` folder
An old attempt at matlab simulation.
TODO(Allen): documentation


# Contribution
By far the easiest way to get up to speed is to read through the `src/ShrimpSimulator.ipynb`. It has detailed notes (although nonlinear) with key equations pulled out for later. However, before jumping into the codebase, it is wise to note the notation used for vector quantities:

### Code Vector Notation:
* Position `r_a2b_c` means the vector from a to b in the c frame.
* Velocity `vel_a2b_c` means the vector of b velocity relative to a velocity in the c frame. For example: `vel_w2b_w` means that since the world is stationary, that's just the b velocity in the world frame.
* Euler angles `eul_a2b` means the euler rotation from a to b, so if you apply it it transforms the a-vector into the b frame. Of course, this does not note the euler angle convention.
* Angular velocity `angvel_a2b_c` means the vector of the b angular velocity relative to a in the c frame. So `angvel_w2b_f` means the angular velocity of the body relative to the world frame (hint, world velocity is zero) in the flyer frame. `w2b_b` is the same as traditional `p,q,r` angular velocities, where the velocities are defined about the body axes of the vehicle.
* Rotation `rot_a2b` means a rotation from a frame to b frame. Applied to an a-frame vector, transforms it into a b frame vector.

### Euler Angles Convention
For euler angles, we are using the `scipy` library of `spatial.transform` with `Rotation`. https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html#scipy.spatial.transform.Rotation

This uses the capital notation `ZXY` to denote intrinsic rotations, which we will use such that the final rotation is about the rotational axis of the vehicle for readibility. This is such that the first rotation is about the Z axis, or corresponds to the yaw.
This matches the MEAM 620 project 1 convention.
Additionally, the euler angles are $\phi, \theta, \psi$, where $\phi$ is yaw.

Please use the `euler2Rotm` and `rotm2Euler` functions for converting between rotation matricies and euler angle vectors.


## Tasks to Do:
* Set up CI
* Write better tests - turn all testing modules into errors
* Add GIF output for shrimpVisualizer
* Hook up inflow velocity
