## ----setup, echo=FALSE---------------------------------------------------
library(knitr)
opts_chunk$set(message=FALSE)

## ----import--------------------------------------------------------------
library(simmer)
library(dplyr)
options(stringsAsFactors = F)

# Load reaction equations (with embedded omic data)
rMech_summary_table <- suppressMessages(readr::read_tsv(system.file("extdata", "reactionEqn_fitting", "rMech_summary_table.tsv", package = "simmer")))
reactionForms <- suppressMessages(readr::read_rds(system.file("extdata", "reactionEqn_fitting", "reactionForms.Rds", package = "simmer")))

kable(head(rMech_summary_table, 8))

## ------------------------------------------------------------------------
# toy run parameters (so everything runs quick)
markov_pars <- list()
markov_pars$sample_freq <- 20 #what fraction of markov samples are reported (this thinning of samples decreases sample autocorrelation)
markov_pars$n_samples <- 200 #how many total markov samples are desired
markov_pars$burn_in <- 0 #how many initial samples should be skipped

## ----evaluate_one_rMech--------------------------------------------------
rxnSummary <- reactionForms[[1]]
kinetically_differing_isoenzymes <- F

# format reaction equations to work with log-concentrations of metabolites and enzymes
rxnEquations <- format_raw_equations(rxnSummary, kinetically_differing_isoenzymes)

# Describe the relevent kinetic parameters

summarize_kinetic_parameters(rxnSummary, rxnEquations, kinetically_differing_isoenzymes)

# Format fluxes, metabolites and enzymes

omic_data <- format_omic_data(kineticPars, all_species, rxnSummary, kinetically_differing_isoenzymes)

# Combine enzyme(s) with non-linear portion of the kinetic form to generate final kinetic form

finalize_reaction_equations(rxnEquations, all_species, kinetically_differing_isoenzymes)

# Formulate priors on non-linear kinetic parameters (not kcat)

kineticParPrior <- build_kinetic_parameter_priors(rxnSummary, kineticPars, omic_data)

# Sampling l(par|X) using MCMC-NNLS

fit_reaction_equation_MCMC_NNLS(markov_pars, kineticPars, kineticParPrior, rxnEquations, all_species, omic_data, kinetically_differing_isoenzymes)

# display posterior samples  
kable(head(markov_par_vals, 10))

