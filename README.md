# Systematic Identification of Meaningful Metabolic Enzyme Regulation (SIMMER)

This repository contains scripts and data that can be used to reproduce the flux inference and kinetic model fitting results of Hackett 2016.

To install this repository as a local package run the following R code:

```{r}
install.packages("devtools")
devtools::install_github("shackett/simmer", build_vingettes = TRUE)
```

Vignettes are available that reproduce the flux inference and reaction equation fitting:

```browseVignettes(package = "simmer")```

## Estimating metabolism-wide fluxes from boundary fluxes using quadratic programming

Integrating boundary fluxes with the yeast metabolic reconstruction in order to estimate metabolic fluxes. This script makes use of the [Gurobi optimizer and R package](http://www.gurobi.com/) which are freely available for academic use.

- Load files describing valid reactions, species (their composition) both from the core SBML model and supplemented manual annotations
- Load files describing boundary conditions, reaction reversibility and auxotrophies
- Setup matrices defining the stoichiometry of each reaction and how reactions will be constrained
- Calculate fluxes using quadratic programming (calculate_QP_fluxes)
- Flux uncertainty can be calculated using python script (python/qp_fba_clust.py)

## R/MCMC_NNLS_example.R

Bayesian approach to fitting reaction equations to measured fluxes.

- Load reaction equations (with embedded omic data)
- Set run parameters (number of samples, sampling frequency and burn-in)
- Fit reaction equations (fit_reaction_equations)
