% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMC_NNLS_functions.R
\name{lik_calc_fittedSD}
\alias{lik_calc_fittedSD}
\title{Find the likelihood of data given proposed parameters}
\usage{
lik_calc_fittedSD(proposed_params, kineticPars, all_species, rxnEquations,
  omic_data, kinetically_differing_isoenzymes)
}
\value{
a numeric log-likelihood
}
\description{
Determine the likelihood of predicted flux as a function of metabolite abundance and
kinetics parameters relative to actual flux
}
