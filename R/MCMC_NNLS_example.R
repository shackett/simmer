options(stringsAsFactors = F)
source("./R/MCMC_NNLS_functions.R")

# Load reaction equations (with embedded omic data)
rMech_summary_table <- read.delim("companionFiles/reactionEqn_fitting/rMech_summary_table.tsv")
reactionForms <- readRDS("companionFiles/reactionEqn_fitting/reactionForms.Rds")

# real run param
#markov_pars <- list()
#markov_pars$sample_freq <- 300 #what fraction of markov samples are reported (this thinning of samples decreases sample autocorrelation)
#markov_pars$n_samples <- 200 #how many total markov samples are desired
#markov_pars$burn_in <- 8000 #how many initial samples should be skipped

# toy run parameters (so everything runs quick)
markov_pars <- list()
markov_pars$sample_freq <- 20 #what fraction of markov samples are reported (this thinning of samples decreases sample autocorrelation)
markov_pars$n_samples <- 200 #how many total markov samples are desired
markov_pars$burn_in <- 0 #how many initial samples should be skipped
 
fitted_reaction_equations <- parallel::mclapply(reactionForms, function(x){
  fit_reaction_equations(x)
}, mc.cores = 4)

# 
rMech_summary_table$maxLogLik <- sapply(fitted_reaction_equations, function(x){max(x[[1]]$logLikelihood)})

library(dplyr)
rMech_summary_table %>%
  group_by(reaction) %>%
  arrange(desc(maxLogLik)) %>% View
