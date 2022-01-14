# Step use-case

This use-case corresponds to the non-parallele Kuramoto-Sivashinsky (KS) equation described in the following article:

> [1] Sipp, D. & Fosas de Pando, M. & Schmid, P.J. (2020). Non-linear model reduction: A comparison between POD-Galerkin and POD-DEIM methods. Computers and Fluids.

where the following statement is made

> This equation is commonly used as a model equation for complex fluid motion as it contains many features in a one-dimensional setting that have equivalents in higher dimensions. It is thus a valuable proxy for investigating analytical and computational techniques and for quantifying the influence of its ingredients on user-specified performance measures.

The spartial discretisation of the PDE leads to a large autonomous (no actuator here) quadratic ODE in limit of stability. Note that the PDE has 1 spatial dimension but is solved within a 2D mesh. This implies in particular that several dof are actually fixed. Routines are provided to get the meaningful parameters.

This benchmark was first made available [on Denis Sipp's Github page](https://github.com/denissipp/CompFluids_SippFosasSchmid_2020/tree/master/KS). The version here is mainly a wrapper to simplify its use for identification/reduction.

## Table of content


## Getting started


```python
import numpy as np
import matplotlib.pyplot as plt
#
import Flocon

# Load the benchmark
ks = Flocon.get_bench('ks')
# Mesh generation
ks.deltax = 0.05 # granularity of the spatial discretisation along the line  
ks.make_mesh()
# compute initial state associated with the 2D-mesh
x0 = ks.make_x0()
```

Note that at this step, the variable `x0` is the one actually used by FreeFem, meaning that is also contains the fixed DOF. To have a graphical representation of a part the system, you may compute the output associated with the (initial) state

```python
[y0 loc_y] = ks.get_y(x0)  
#
plt.plot(loc_y, y0)
plt.show()
```

Based on the initial state, you may also launch a simulation which returns the time, the energy of the system and the output. This is done as follows,
```python
[t, energy, y]= ks.simulate(x0)
```

## Acknowledgement  

Special thanks to Denis Sipp for the use-case.
