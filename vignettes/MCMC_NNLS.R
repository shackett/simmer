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

