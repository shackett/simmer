# Systematic Identification of Meaningful Metabolic Enzyme Regulation (SIMMER)

This repository contains scripts and data that can be used to reproduce the flux inference and kinetic model fitting results of Hackett 2016.

## R/QP_FBA_example.R

Integrating boundary fluxes with the yeast metabolic reconstruction in order to estimate metabolic fluxes. THis script makes use of the [Gurobi optimizer](http://www.gurobi.com/).

- Load files describing valid reactions, species (their composition) both from the core SBML model and supplemented manual annotations
- Load files describing boundary conditions, reaction reversibility and auxotrophies
- Setup matrices defining the stoichiometry of each reaction and how reactions will be constrained
- Calculate fluxes using quadratic programming
