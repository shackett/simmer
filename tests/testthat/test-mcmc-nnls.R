test_that("fit_reaction_equations applies the MCMC-NNLS algorithm to one reaction equation", {

#rMech_summary_table <- read.delim("companionFiles/reactionEqn_fitting/rMech_summary_table.tsv")
#reactionForms <- readRDS("companionFiles/reactionEqn_fitting/reactionForms.Rds")

library(simmer)

data(rMech_summary_table)
data(reactionForms)
expect_true(all(c("rMech_summary_table", "reactionForms") %in% ls()))
expect_equal(colnames(rMech_summary_table), c("reaction", "rMech", "modelType"))
expect_equal(rMech_summary_table$rMech, names(reactionForms))

markov_pars <- list()
markov_pars$sample_freq <- 20 #what fraction of markov samples are reported (this thinning of samples decreases sample autocorrelation)
markov_pars$n_samples <- 200 #how many total markov samples are desired
markov_pars$burn_in <- 0 #how many initial samples should be skipped

fitted_reaction_equations <- fit_reaction_equations(reactionForms[[1]])

expect_equal(names(fitted_reaction_equations), names(reactionForms[1]))
expect_equal(names(fitted_reaction_equations[[1]]), c("kineticPars", "all_species", "kineticParPrior", "markovChain",
                                                      "logLikelihood", "occupancyEq", "metabolites", "enzymes", "flux",
                                                      "specSD", "specCorr"))
expect_equal(nrow(fitted_reaction_equations[[1]]$markovChain), markov_pars$n_samples)
})
