# PDEfind

A Julia implementation of the *PDE functional identification of nonlinear dynamics* (PDEfind) algorithm, first introduced in \\
Samuel H. Rudy *et al.*,  
Data-driven discovery of partial differential equations.
*Sci. Adv.* **3**, e1602614 (2017).  
[DOI: 10.1126/sciadv.1602614](https://doi.org/10.1126/sciadv.1602614)
 
## Outline (and todo)

### Data collection
- [x] Canonical PDEs (only 1d heat eq.)
- [ ] Add noise
- [ ] Visualisation framework

### Build library of data and derivatives
- [ ] Available functions selection
- [x] Derivative degree selection
- [ ] Derivation method selection
- [x] Preprocess library matrix columns to unit variance
- [x] Build library of data and derivatives

### Solve sparse regression
- [x] Select tolerance, regularization factor and number of iterations
- [x] Run Sequential Threshold Ridge regression
- [x] Print selected PDE and coefficients


[![Build Status](https://github.com/dknatan/PDEfind.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/dknatan/PDEfind.jl/actions/workflows/CI.yml?query=branch%3Amain)
