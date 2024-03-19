# PDEfind

A Julia implementation of the *PDE functional identification of nonlinear dynamics* (PDEfind) algorithm, first introduced in  Samuel H. Rudy *et al.*,  
Data-driven discovery of partial differential equations. *Sci. Adv.* **3**, e1602614 (2017). [DOI: 10.1126/sciadv.1602614](https://doi.org/10.1126/sciadv.1602614).
This version only works for 1 component PDEs in 1 spatial dimension.

## Outline

### Mock data preparation
Mock data prepared with the following packages: `OrdinaryDiffEq, ModelingToolkit, MethodOfLines, DomainSets`. 

### Discretization
The grid along with finite difference derivative operators, that serve as the inputs for the function libraries.

### Build function library
For now only the polynomial basis is implemented. It constructs a library of polynomials of a specified maximum derivative and degree.

### Solve sparse regression
Both algorithms from Rudy et al. are implemented and can be called directly. 

[![Build Status](https://github.com/dknatan/PDEfind.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/dknatan/PDEfind.jl/actions/workflows/CI.yml?query=branch%3Amain)
