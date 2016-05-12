populate_reactionEqns <- function(){

# Not meant for users
# Setup a subset of reaction forms for package
# including rMM reaction forms, significant complex regulation (2 regulators or cooperativity) and
# either significant single regulors (or best regulators for complex regulation)
# also a couple of examples showing that the method can be used to test isoenzyme-specific kinetics (regulator or metabolite affinities)

rxn_forms <- all_rxn_fits %>% arrange(desc(modelType))

rxn_forms <- rbind(rxn_forms,
  data.frame(reaction = "r_0042",
                        rMech = c("r_0042_Y_F_inhibition_isoenzymeSpecific", "r_0042_E4P_enzyme_specific_affinity_test2"),
                        modelType = c("example: testing isoenzyme-specific regulation", "example: testing isoenzyme-specific substrate affinities")))

write.table(rxn_forms, file = "companionFiles/reactionEqn_fitting/rMech_summary_table.tsv", sep = "\t", row.names = F, quote = F)

load("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale/paramOptim.Rdata")

reactionForms <- rxnList_form[rxn_forms$rMech]

saveRDS(reactionForms, file = "companionFiles/reactionEqn_fitting/reactionForms.Rds")

}


source("R/MCMC_NNLS_functions.R")

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
markov_pars$sample_freq <- 5 #what fraction of markov samples are reported (this thinning of samples decreases sample autocorrelation)
markov_pars$n_samples <- 10 #how many total markov samples are desired
markov_pars$burn_in <- 0 #how many initial samples should be skipped
  
rxnSummary <- reactionForms[[1]]


fit_reactionEqn <- function(rxnSummary)

  ##### Generate non-linear part(s) of kinetic form ####
  
  # this method allows for isoenzymes with different kinetics (none of the published work used this feature due to the extra
  # degrees of freedom it uses
  kinetically_differing_isoenzymes <- any(names(rxnSummary$rxnForm) %in% rownames(rxnSummary$enzymeComplexes))
  # if isoenzymes differ w.r.t. kinetics or regulation then their occupancy equation are stored as different elements of a list
  # and enzyme concentrations need to be paired with these seperate equations
  
  # format reaction equations to work with log-concentrations of metabolites and enzymes
  rxnEquations <- format_raw_equations(rxnSummary, kinetically_differing_isoenzymes)

  # Describe the relevent kinetic parameters
  
  summarize_kinetic_parameters(rxnSummary, kinetically_differing_isoenzymes)
  
  # Format fluxes, metabolites and enzymes
  
  omic_data <- format_omic_data(kineticPars, all_species, kinetically_differing_isoenzymes)
  
  # Combine enzyme(s) with non-linear portion of the kinetic form to generate final kinetic form
  
  finalize_reaction_equations(rxnEquations, all_species, kinetically_differing_isoenzymes)
  
  # Formulate priors on non-linear kinetic parameters (not kcat)
  
  n_c <- nrow(omic_data$flux)
  
  kineticParPrior <- build_kinetic_parameter_priors(kineticPars)
  
  #### Sampling l(par|X) using MCMC-NNLS ####
  
  fit_reaction_equation_MCMC_NNLS(markov_pars, omic_data, kineticParPrior)
  
  

  ## Save MC information
  run_summary[[names(rxnList_form)[rxN]]]$kineticPars <- kineticPars
  run_summary[[names(rxnList_form)[rxN]]]$markovChain <- markov_par_vals
  run_summary[[names(rxnList_form)[rxN]]]$likelihood <- lik_track
  ##  Save flux and rxn information
  run_summary[[names(rxnList_form)[rxN]]]$metabolites <- met_abund
  run_summary[[names(rxnList_form)[rxN]]]$enzymes <- enzyme_abund
  run_summary[[names(rxnList_form)[rxN]]]$flux <- flux
  run_summary[[names(rxnList_form)[rxN]]]$occupancyEq <- rxnEquations
  run_summary[[names(rxnList_form)[rxN]]]$all_species <- all_species
  run_summary[[names(rxnList_form)[rxN]]]$rxnSummary <- rxnSummary
  run_summary[[names(rxnList_form)[rxN]]]$specSD <- species_SD
  run_summary[[names(rxnList_form)[rxN]]]$specCorr <- species_corr
  ## Save priors
  run_summary[[names(rxnList_form)[rxN]]]$kineticParPrior <- data.frame(kineticPars, kineticParPrior)
  
  


