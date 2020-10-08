# Step use-case

This use-case corresponds to the benchmark described in the following article

> [1] Sipp, D. & Schmid, P. J. (2016). Linear Closed-Loop Control of Fluid Instabilities and Noise-Induced Perturbations: A Review of Approaches and Tools. Applied Mechanics Reviews.

This benchmark was first made available [on Denis Sipp's Github page](https://github.com/denissipp/AMR_Sipp_Schmid_2016). The version here is mainly a wrapper than simplifies its use for identification and control design and allows for a simple modification of the input/output configuration.

## Table of content

* [Introduction](#introduction)
  * [Requirements](#requirements)
* [Getting started](#getting-started)
  * [Initialisation](#initialisation)
  * [Base-flow](#base-flow-computation)
  * [Open-loop simulation](#open-loop-simulation)
  * [Closed-loop simulation](#closed-loop-simulation)
* [Modifying the configuration](#modifying-the-configuration)  
* [Acknowledgement](#acknowledgement)

## Introduction

The default configuration is the following:
![Default configuration](./static/step_default.png)

## Requirements

* [FreeFem++](https://freefem.org/) for the core of the simulator
* numpy (>= 1.18) for input/output signals
* matplotlib


## Getting started


### Initialisation

The benchmark can be loaded as follows,

```python
import numpy as np
import matplotlib.pyplot as plt
#
import Flocon
# Load the step benchmark
step = Flocon.get_bench('step')
# By default, the configuration is the same as in [1]
step.plot_config()
```
This creates a (sqlite) database in which simulations results are stored.


### Base-flow computation

Simulations are performed around an equilibrium called base-flow which  must be computed as follows:

```python
# Base-flow computation
step.Re   = 100     # change Reynolds number
step.mesh = 'test'  # change active mesh  
step.compute_baseflow()
```
The results is directly stored in the database. The base-flow depends both on the Reynolds number and on the mesh. Changing one of those parameters requires to computer another base-flow.

*About convergence of the base-flow:* Note that the computation of the base-flow is done with a Newton iteration which may not converge, especially for large Reynolds number. To alleviate this issue, an automated procedure has been implemented: attempt to solve for the current Reynolds number by restarting from the closest existing solution (stored in DB). If it converges, store the solution. Otherwise, decrease the Reynolds number and repeat.

However, this trial and error scheme may take a while and you can also manually increase the Reynolds number.

*About the meshes:* Three increasingly refined mesh can be considered: `'test'`, `'coarse'` and `'fine'`. The test mesh is inaccurate and should be used, as its name suggests, only to play around with the benchmark. For instance, it may help understand what the various routines are expecting and returning. On the contrary, the other two meshes are way more refined and induce a significant computational burden for simulation.

### Open-loop simulation

Once the base-flow has been computed, open-loop simulations can be performed. Those, in addition to the Reynolds number and mesh, depends on several additional parameters listed below.

|   Name              | Description                                            |  Default value |
|--------------------:|:-------------------------------------------------------|:--------------:|
| `step.dt`           | time step                                              |  `.002`        |
| `step.N`            | number of time samples                                 |  `250000`      |
| `step.NL`           | flag indicating whether the non-linearity is activated |  `False`       |
| `step.noises`       | list of perturbations                                  |  as in [1]     |
| `step.actuators`    | list of control inputs                                 |  as in [1]     |
| `step.performances` | list of performances outputs                           |  as in [1]     |
| `step.sensors`      | list of sensors                                        |  as in [1]     |

*Remark: `N` can be modified indirectly by setting the final time `step.tf`. The time step remains untouched.*

Note that in addition to the physical configuration of the inputs and outputs (see the associated sections below), an open-loop simulation depends on the exact signal that is fed to the system.
An open-loop simulation is therefore stored in the database as the parameters in the table above, the matrices containing the input and output signals and a name (to make it intelligible). The latter is used as identifier to retrieve the data from the database.

By default, perturbations are handled within Flocon so that the user need only to provide the control input. For instance, to simulate the open-loop with no command input (but with default noise):

```python
u       = np.zeros([step.N,1])              # control signal
name    = 'OL1'                             # name of the simulation
out     = step.simulate_openloop(name, u)   # run the simulation
t       = out[:,0]                          # time is the first output
energy  = out[:,1]                          # then comes the flow energy
z       = out[:,2]                          # then comes the performances outputs
y       = out[:,3]                          # and the measurement outputs
# plotting the signals
plt.plot(t, z, t, y)
plt.show()
```
If you attempt to make a simulation which name already exists in the database, the simulation outputs are returned. To get both input (noise included) and output data, use `step.get_openloop_simulation(name)`.

*Remark: For an experiment to be reproducible, you can fix the seed of the random noise generator by providing it `step.simulate_openloop('OL', u, seed=<user_seed>)`.*

### Closed-loop simulation

Assuming that you have a controller `K` (a dictionary) described by a *discrete* linear state-space representation `A`, `B`, `C` and `D`, then you can launch a closed-loop simulation as follows:

```python
# K is a dictionary with fields A, B, C, D
name    = 'CL1'                             # name of the closed-loop
out     = step.simulate_closedloop(name, K) # run the simulation
t       = out[:,0]                          # time is the first output
energy  = out[:,1]                          # then comes the flow energy
zcl     = out[:,2]                          # then comes the performances outputs
ycl     = out[:,3]                          # and the measurement outputs
# plotting the signals
plt.plot(t, zcl, t, ycl)
plt.show()
```

By default, the control-law is placed in a positive feedback. To add a reference and feed the control-law with a tracking error, it must be provided when launching the simulation as `out = step.simulate_closedloop(name, K, ref = r)` where `r` is a numpy array containing the reference.

As with the open-loop simulation, the input/output signals and the control-law data can be recovered as a dictionary with `data = step.get_closedloop_simulation(name)`. If a reference was added, then it appears in the last column of the field `in`.

## Modifying the configuration

WIP

## Acknowledgement  

Special thanks to Denis Sipp and Colin Leclercq for their help with the integration of this use-case.
