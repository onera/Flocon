# Step use-case

This use-case corresponds to the benchmark described in the following article

> [1] Hervé, A., Sipp, D., Schmid, P. J., & Samuelides, M. (2012). A physics-based approach to flow control using system identification. Journal of Fluid Mechanics, 702, 26-58.


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
# step.Re   = 100
# step.mesh ='coarse' # other possibility is 'fine'
# Base-flow computation
step.compute_baseflow()
```
The results is directly stored in the database. The base-flow depends both on the Reynolds number and on the mesh. Changing one of those parameters requires to computer another base-flow.


*About convergence of the base-flow:* Note that the computation of the base-flow is done with a Newton iteration which may not converge, especially for large Reynolds number. To alleviate this issue, an automated procedure has been implemented: attempt to solve for the current Reynolds number by restarting from the closest existing solution (stored in DB). If it converges, store the solution. Otherwise, decrease the Reynolds number and repeat.

However, this trial and error scheme may take a while and you can also manually increase the Reynolds number.

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
out     = step.simulate_openloop('OL', u)   # name of the simulation and command
t       = out[:,0]                          # time is the first output
z       = out[:,1]                          # then comes the performances outputs
y       = out[:,2]                          # and the measurement outputs
# plotting the signals
plt.plot(t, z, t, y)
plt.show()
```
If you attempt to make a simulation which name already exists in the database, the simulation outputs are returned. To get both input and output data, use `step.get_openloop_simulation(name)`.

*Remark: For an experiment to be reproducible, you can fix the seed of the random noise generator by providing it `step.simulate_openloop('OL',u, seed=<user_seed>)`.*

### Closed-loop simulation

Assuming that you have a controller `K` (a dictionary) described by a *discrete* linear state-space representation `A`, `B`, `C` and `D`, then you can launch a closed-loop simulation as follows:

```python
# K is a dictionary with fields A, B, C, D
out     = step.simulate_closedloop('CL', K) # name of the simulation and control law
t       = out[:,0]                          # time is the first output
zcl     = out[:,1]                          # then comes the performances outputs
ycl     = out[:,2]                          # and the measurement outputs
# plotting the signals
plt.plot(t, zcl, t, ycl)
plt.show()
```



## Modifying the configuration

In order to perform some actuator or sensor placement studies, it is possible to add/modify the inputs and outputs for this configuration.


## Acknowledgement  

Special thanks to Denis Sipp and Colin Leclercq for their precious help for the integration of this use-case.
