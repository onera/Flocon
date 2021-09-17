# Step use-case

This use-case corresponds to the non-parallele Kuramoto-Sivashinsky (KS) equation described in the following article:

> [1] Sipp, D. & Fosas de Pando, M. & Schmid, P.J. (2020). Non-linear model reduction: A comparison between POD-Galerkin and POD-DEIM methods. Computers and Fluids.

where the following statement is made

> This equation is commonly used as a model equation for complex fluid motion as it contains many features in a one-dimensional setting that have equivalents in higher dimensions. It is thus a valuable proxy for investigating analytical and computational techniques and for quantifying the influence of its ingredients on user-specified performance measures.

The spartial discretisation of the PDE leads to an autonomous (no actuator here) quadratic ODE of dimension 16 000. The linearised model is unstable while the overall quadratic model is stable.

This benchmark was first made available [on Denis Sipp's Github page](https://github.com/denissipp/CompFluids_SippFosasSchmid_2020/tree/master/KS). The version here is mainly a wrapper than simplifies its use for identification/reduction.

## Table of content



## Acknowledgement  

Special thanks to Denis Sipp for the use-case.
