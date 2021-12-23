# PhD_Codes

This repository contains codes developed during my PhD at University of Groningen.
- [OpenFOAM](/OpenFOAM) contains applications (solvers, utilities) and src (fvOptions) compiled with version (v2106).
- [scripts](/scripts) contains Python and MATLAB scripts for pre- and post- processing. Also contains TensorFlow modular code for data-driven project for UQ analysis.

## [OpenFOAM](/OpenFOAM)

### Introduction
Using an exiting deterministec solver ```pimpleFoam```, a stoachastic solver based on Intrusive Polynomial Chaos (IPC) was developed in OpenFOAM called as [```gPCModelFormSimpleFoam```](/OpenFOAM/p285464-v2012/applications/solvers/). More details on this solver can be found in the following articles:

- Quantification and propagation of model-form uncertainties in RANS turbulence modeling via intrusive polynomial chaos. _International Journal for Uncertainty Quantification_. (under review)
- Intrusive Polynomial Chaos for CFD Using OpenFOAM. Computational Science, _ICCS 2020, Lecture Notes in Computer Science. Springer, Cham._

### Stochatic Simulation Setup
