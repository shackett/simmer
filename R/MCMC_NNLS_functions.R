fit_reaction_equations <- function(rxnSummary){
  
  # This function takes one reaction equation and associated multi-omic data and evaluate how well the
  # reaction equation can predict provided flux based on the MCMC-NNLS algorithm.
  
  # this method allows for isoenzymes with different kinetics (none of the published work used this feature due to the extra
  # degrees of freedom it uses
  kinetically_differing_isoenzymes <- any(names(rxnSummary$rxnForm) %in% rownames(rxnSummary$enzymeComplexes))
  # if isoenzymes differ w.r.t. kinetics or regulation then their occupancy equation are stored as different elements of a list
  # and enzyme concentrations need to be paired with these seperate equations
  
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
  
  #### Sampling l(par|X) using MCMC-NNLS ####
  
  fit_reaction_equation_MCMC_NNLS(markov_pars, kineticPars, kineticParPrior, rxnEquations, all_species, omic_data, kinetically_differing_isoenzymes)
  
  run_summary <- list()
  ## Save MC information
  run_summary[[rxnSummary$listEntry]]$kineticPars <- kineticPars
  run_summary[[rxnSummary$listEntry]]$all_species <- all_species
  run_summary[[rxnSummary$listEntry]]$kineticParPrior <- data.frame(kineticPars, kineticParPrior)
  run_summary[[rxnSummary$listEntry]]$markovChain <- markov_par_vals
  run_summary[[rxnSummary$listEntry]]$logLikelihood <- lik_track
  ##  Save flux and rxn information
  run_summary[[rxnSummary$listEntry]]$occupancyEq <- rxnEquations
  run_summary[[rxnSummary$listEntry]]$metabolites <- omic_data$met_abund
  run_summary[[rxnSummary$listEntry]]$enzymes <- omic_data$enzyme_abund
  run_summary[[rxnSummary$listEntry]]$flux <- omic_data$flux
  run_summary[[rxnSummary$listEntry]]$specSD <- omic_data$species_SD
  run_summary[[rxnSummary$listEntry]]$specCorr <- omic_data$species_corr
  return(run_summary)
}

#### Setting up reaction equations and omic data
 
format_raw_equations <- function(rxnSummary, kinetically_differing_isoenzymes){
  
  # reformat reactionEquations to work with log-space data
  
  rxnEquations <- list()
  
  if(kinetically_differing_isoenzymes){
    for(isoenzyme in names(rxnSummary$rxnForm)){
      
      rxnEquations$occupancyEq_list[[isoenzyme]] <- paste(deparse(as.list(rxnSummary$rxnForm[names(rxnSummary$rxnForm) == isoenzyme])[[1]][[2]]), collapse = "") # a parametric form relating metabolites and constants to fraction of maximal activity
      rxnEquations$occupancyEq_list[[isoenzyme]] <- gsub('[ ]+', ' ', rxnEquations$occupancyEq_list[[isoenzyme]])
      
      rxnEquations$l_occupancyEq_list[[isoenzyme]] <- rxnEquations$occupancyEq_list[[isoenzyme]]
      rxnEquations$l_occupancyEq_list[[isoenzyme]] <- gsub('([^_])(t_)', '\\12^\\2', rxnEquations$l_occupancyEq_list[[isoenzyme]])
      rxnEquations$l_occupancyEq_list[[isoenzyme]] <- gsub('([0-9]+)\\^(t_[0-9]{4})\\^([0-9]+)', '(\\1\\^\\2)\\^\\3', rxnEquations$l_occupancyEq_list[[isoenzyme]]) # correct 2^X^2 -> (2^X)^2
      
      # create an expression for each isoenzyme
      rxnEquations[["l_occupancyExpression"]][[isoenzyme]] <- parse(text = sub(' \\+ 0$', '', sub('^I', '', rxnEquations$l_occupancyEq_list[[isoenzyme]]))) 
    }
  }else{
    # the typical case - single enzyme or kinetically equivalent isoenzymes
    rxnEquations[["occupancyEq_list"]] <- paste(deparse(as.list(rxnSummary$rxnForm)[[2]]), collapse = "") # a parametric form relating metabolites and constants to fraction of maximal activity
    rxnEquations[["occupancyEq_list"]] <- gsub('[ ]+', ' ', rxnEquations[["occupancyEq_list"]])
    
    rxnEquations[["l_occupancyEq_list"]] <- rxnEquations[["occupancyEq_list"]]
    rxnEquations[["l_occupancyEq_list"]] <- gsub('([^_])(t_)', '\\12^\\2', rxnEquations[["l_occupancyEq_list"]])
    rxnEquations[["l_occupancyEq_list"]] <- gsub('([0-9]+)\\^(t_[0-9]{4})\\^([0-9]+)', '(\\1\\^\\2)\\^\\3', rxnEquations[["l_occupancyEq_list"]]) # correct 2^X^2 -> (2^X)^2
    
    # create an expression for each isoenzyme
    rxnEquations[["l_occupancyExpression"]] <- parse(text = sub(' \\+ 0$', '', sub('^I', '', rxnEquations$l_occupancyEq_list))) 
  }
  
  return(rxnEquations)
  
  }

summarize_kinetic_parameters <- function(rxnSummary, rxnEquations, kinetically_differing_isoenzymes){
  
  # based on metabolites and enzyme data provided, determine the kinetic constants that are present in the reaction equations
  
  if(kinetically_differing_isoenzymes & length(grep('.[0-9]$', names(rxnSummary$rxnForm))) != 0){
    # A single enzyme is applied to multiple expressions (ex. partial constitutive activity)
    kineticPars <- data.frame(rel_spec = c(names(rxnSummary$rxnForm), rxnSummary$rxnFormData$SubstrateID), 
                              SpeciesType = c(rep("Enzyme", times = length(rxnSummary$rxnForm)), rep("Metabolite", times = nrow(rxnSummary$rxnFormData))), modelName = NA, commonName = NA, formulaName = NA, measured = NA)
  }else{
    # One enzyme, one Kcat
    kineticPars <- data.frame(rel_spec = c(rownames(rxnSummary$enzymeComplexes), rxnSummary$rxnFormData$SubstrateID), 
                              SpeciesType = c(rep("Enzyme", times = nrow(rxnSummary$enzymeComplexes)), rep("Metabolite", times = nrow(rxnSummary$rxnFormData))), modelName = NA, commonName = NA, formulaName = NA, measured = NA)
  }
  kineticPars$formulaName[kineticPars$SpeciesType == "Enzyme"] <- paste(paste("E", rxnSummary$rxnID, sep = "_"), kineticPars$rel_spec[kineticPars$SpeciesType == "Enzyme"], sep = "_")
  kineticPars$modelName[kineticPars$SpeciesType == "Metabolite"] <- kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"]
  kineticPars$commonName[kineticPars$SpeciesType == "Metabolite"] <- unname(sapply(kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"], function(x){rxnSummary$metNames[names(rxnSummary$metNames) == x]}))
  kineticPars$commonName[kineticPars$SpeciesType == "Enzyme"] <- kineticPars$rel_spec[kineticPars$SpeciesType == "Enzyme"]
  kineticPars$formulaName[kineticPars$SpeciesType == "Metabolite"] <-  rxnSummary$rxnFormData$affinityParameter
  
  all_species <- kineticPars[sapply(kineticPars$formulaName, function(ele_used){length(grep(ele_used, rxnEquations[["occupancyEq_list"]])) != 0}) | kineticPars$SpeciesType == "Enzyme",]
  kineticPars <- kineticPars[sapply(kineticPars$formulaName, function(ele_used){length(grep(ele_used, rxnEquations[["occupancyEq_list"]])) != 0}),] #remove species which dont appear in the reaction equation
  
  kineticPars <- rbind(kineticPars, c("keq", "keq", NA, NA, paste("Keq", rxnSummary$rxnID, sep = ""), NA))
  if(any(rxnSummary$rxnFormData$Hill == 0)){ # an unspecified hill coefficient was found - introduce hill parameter
    kineticPars <- rbind(kineticPars, c(rxnSummary$rxnFormData$SubstrateID[rxnSummary$rxnFormData$Hill == 0], "hillCoefficient", NA, NA, sub('^K', 'KH', rxnSummary$rxnFormData$affinityParameter[rxnSummary$rxnFormData$Hill == 0]), NA))
  }
  
  kineticPars <<- kineticPars
  all_species <<- all_species
  
}

format_omic_data <- function(kineticPars, all_species, rxnSummary, kinetically_differing_isoenzymes){
  
  ### Create a matrix containing the metabolites and enzymes
  if(kinetically_differing_isoenzymes & length(grep('.[0-9]$', names(rxnSummary$rxnForm))) != 0){
    # A single enzyme is applied to multiple expressions (ex. partial constitutive activity)
    enzyme_abund <- t(rxnSummary$enzymeComplexes)[,data.table::chmatch(gsub('.[0-9]$', '', names(rxnSummary$rxnForm)), rownames(rxnSummary$enzymeComplexes))]
    colnames(enzyme_abund) <- names(rxnSummary$rxnForm)
    
    species_SD <- rxnSummary$all_species_SD
    remapping_indices <- data.table::chmatch(gsub('.[0-9]$', '', names(rxnSummary$rxnForm)), colnames(species_SD)[colnames(species_SD) %in% rownames(rxnSummary$enzymeComplexes)])
    remapping_table <- data.frame(remap = names(rxnSummary$rxnForm), enzyme = names(remapping_indices))
    renamed_enzyme_SD <- species_SD[,colnames(species_SD) %in% rownames(rxnSummary$enzymeComplexes)][,remapping_indices]
    colnames(renamed_enzyme_SD) <- names(rxnSummary$rxnForm)
    species_SD <- cbind(species_SD[,!(colnames(species_SD) %in% rownames(rxnSummary$enzymeComplexes))], renamed_enzyme_SD)
    
    species_corr <- matrix(0, nrow = ncol(species_SD), ncol = ncol(species_SD))
    rownames(species_corr) <- colnames(species_corr) <- colnames(species_SD)
    species_corr[!(colnames(species_SD) %in% names(rxnSummary$rxnForm)),!(colnames(species_SD) %in% names(rxnSummary$rxnForm))] <- rxnSummary$all_species_corr[!(rownames(rxnSummary$all_species_corr) %in% rownames(rxnSummary$enzymeComplexes)),!(rownames(rxnSummary$all_species_corr) %in% rownames(rxnSummary$enzymeComplexes))]
    
    remapped_enzyme_corr <- species_corr[(colnames(species_SD) %in% names(rxnSummary$rxnForm)),(colnames(species_SD) %in% names(rxnSummary$rxnForm))]
    for(i in 1:nrow(remapped_enzyme_corr)){
      for(j in 1:ncol(remapped_enzyme_corr)){
        remapped_enzyme_corr[i,j] <- ifelse(remapping_table$enzyme[remapping_table$remap == colnames(remapped_enzyme_corr)[j]] == remapping_table$enzyme[remapping_table$remap == rownames(remapped_enzyme_corr)[i]], 1, 0)
      }
    }
    
    species_corr[rownames(species_corr) %in% rownames(remapped_enzyme_corr), colnames(species_corr) %in% colnames(remapped_enzyme_corr)] <- remapped_enzyme_corr    
    
  }else{
    enzyme_abund <- t(rxnSummary$enzymeComplexes); colnames(enzyme_abund) <- all_species$rel_spec[all_species$SpeciesType == "Enzyme"]
    species_SD <- rxnSummary$all_species_SD
    species_corr <- rxnSummary$all_species_corr
    
  }
  
  met_abund <- rxnSummary$rxnMet
  met_abund <- met_abund[,colnames(met_abund) %in% kineticPars$rel_spec]
  
  if(length(kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"]) <= 1){
    met_abund <- data.frame(met_abund)
    colnames(met_abund) <- kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"]
    kineticPars$measured[kineticPars$SpeciesType == "Metabolite"] <- !all(is.na(met_abund))
  }else{
    kineticPars$measured[kineticPars$SpeciesType == "Metabolite"] <- unname(sapply(kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"], function(x){(apply(is.na(met_abund), 2, sum) == 0)[names((apply(is.na(met_abund), 2, sum) == 0)) == x]}))
  }
  kineticPars$measured <- as.logical(kineticPars$measured)
  
  met_abund[,!kineticPars$measured[data.table::chmatch(colnames(met_abund), kineticPars$modelName)]] <- 0 # set missing data (unascertained) to invariant across conditions
  
  flux <- rxnSummary$flux/median(abs(rxnSummary$flux$standardQP[rxnSummary$flux$standardQP != 0])) #flux, scaled to a prettier range
  
  # If FVA min flux is greater than max flux, switch them (and print a warning).
  
  if(sum(!(flux$FVAmax >= flux$FVAmin)) != 0){
    print(paste("maximum flux is less than minimum flux for", sum(!(flux$FVAmax >= flux$FVAmin)), "conditions"))
  }
  flux[!(flux$FVAmax >= flux$FVAmin),c('FVAmin', 'FVAmax')] <- flux[!(flux$FVAmax >= flux$FVAmin),c('FVAmax', 'FVAmin')]
  
  # If bounds are exactly equal, then introduce a minute spread so a range can be calculated ###
  
  flux$FVAmax[flux$FVAmax == flux$FVAmin] <- flux$FVAmax[flux$FVAmax == flux$FVAmin] + abs(flux$FVAmax[flux$FVAmax == flux$FVAmin]*10^-4)
  
  # Assign some objects to reaction_equations environment
  kineticPars <<- kineticPars
  
  omic_data <- list()
  omic_data$enzyme_abund = enzyme_abund
  omic_data$met_abund = met_abund
  omic_data$species_SD = species_SD
  omic_data$species_corr = species_corr
  omic_data$flux = flux
  return(omic_data)
  
}

finalize_reaction_equations <- function(rxnEquations, all_species, kinetically_differing_isoenzymes){
  
  # expression combining the log-occupancy equation and scaling of enzyme abundance by activity
  # determined w.r.t. metabolites/enzymes in log-space and linear-space
  # if there are multiple isoenzymes which differ kinetically then their expressions are generated seperately and then concatenated
  # a single equation is returned
  
  if(kinetically_differing_isoenzymes){
    # if isoenzymes differ then their occupancy equations are not distributed e.g. Vmax1[O1] + Vmax2[O2] rather than (Vmax1 + Vmax2)*O
    KcatEs_log <- mapply(function(E, Kcat){paste(Kcat, E, sep = " * ")}, E = sapply(all_species$rel_spec[all_species$SpeciesType == "Enzyme"], function(x){paste("2^", x, sep = "")}), Kcat = all_species$formulaName[all_species$SpeciesType == "Enzyme"])
    KcatEs_linear <- mapply(function(E, Kcat){paste(Kcat, E, sep = " * ")}, E = all_species$rel_spec[all_species$SpeciesType == "Enzyme"], Kcat = all_species$formulaName[all_species$SpeciesType == "Enzyme"])
    
    if(!all(names(rxnEquations$occupancyEq_list) %in% c(names(KcatEs_log), "other"))){
      stop(paste("isoenzyme name does not match available complexes for reaction", names(rxnList_form)[rxN],"\nCheck the field \"enzymeInvolved\" in manual_ComplexRegulation.tsv"))
    }
    
    # a complex is matched to an enzyme either by being directly specified or otherwise being placed with "other" enzymes 
    occEqtn_complex_match <- data.frame(complex = names(KcatEs_log), occEqtn = NA)
    occEqtn_complex_match$occEqtn[occEqtn_complex_match$complex %in% names(rxnEquations$occupancyEq_list)] <- occEqtn_complex_match$complex[occEqtn_complex_match$complex %in% names(rxnEquations$occupancyEq_list)]
    occEqtn_complex_match$occEqtn[is.na(occEqtn_complex_match$occEqtn)] <- "other"
    
    KcatExpressions_linear <- NULL
    KcatExpressions_log <- NULL
    for(isoenzyme in unique(occEqtn_complex_match$occEqtn)){
      # log
      KcatExpression <- paste('(', paste(KcatEs_log[names(KcatEs_log) %in% occEqtn_complex_match$complex[occEqtn_complex_match$occEqtn == isoenzyme]], collapse = " + "), ')', sep = "")
      isoenzymeO <- rxnEquations$l_occupancyEq_list[names(rxnEquations$occupancyEq_list) == isoenzyme]
      isoenzymeO <- sub('^I', paste(KcatExpression, '*', sep = ''), isoenzymeO)
      isoenzymeO <- sub(' \\+ 0$', '', isoenzymeO)
      
      KcatExpressions_log <- c(KcatExpressions_log, isoenzymeO)
      
      # linear
      KcatExpression <- paste('(', paste(KcatEs_linear[names(KcatEs_linear) %in% occEqtn_complex_match$complex[occEqtn_complex_match$occEqtn == isoenzyme]], collapse = " + "), ')', sep = "")
      isoenzymeO <- rxnEquations$occupancyEq_list[names(rxnEquations$occupancyEq_list) == isoenzyme]
      isoenzymeO <- sub('^I', paste(KcatExpression, '*', sep = ''), isoenzymeO)
      isoenzymeO <- sub(' \\+ 0$', '', isoenzymeO)
      
      KcatExpressions_linear <- c(KcatExpressions_linear, isoenzymeO)
    }
    
    rxnEquations[["full_kinetic_form_list"]] <- paste(unlist(KcatExpressions_log), collapse = " + ")
    rxnEquations[["elasticity_calc"]] <- paste(unlist(KcatExpressions_linear), collapse = " + ")
    
  }else{
    
    # (Vmax1 + Vmax2)*O
    
    KcatEs <- mapply(function(E, Kcat){paste(Kcat, E, sep = " * ")}, E = sapply(all_species$rel_spec[all_species$SpeciesType == "Enzyme"], function(x){paste("2^", x, sep = "")}), Kcat = all_species$formulaName[all_species$SpeciesType == "Enzyme"])
    KcatExpression <- paste('(', paste(KcatEs, collapse = " + "), ')', sep = "")
    Kcatpaste <- paste('I(', KcatExpression, ' * ', sep = "")
    
    rxnEquations[["full_kinetic_form_list"]] <- rxnEquations[["l_occupancyEq_list"]]
    rxnEquations[["full_kinetic_form_list"]] <- gsub('(I\\()', Kcatpaste, rxnEquations[["full_kinetic_form_list"]])
    rxnEquations[["full_kinetic_form_list"]] <- sub('I\\(', '\\(', rxnEquations[["full_kinetic_form_list"]])
    
    # save the same expression in linear space, so that it can be used later on to look at elasticitity
    
    KcatEs <- mapply(function(E, Kcat){paste(Kcat, E, sep = " * ")}, E = all_species$rel_spec[all_species$SpeciesType == "Enzyme"], Kcat = all_species$formulaName[all_species$SpeciesType == "Enzyme"])
    KcatExpression <- paste('(', paste(KcatEs, collapse = " + "), ')', sep = "")
    Kcatpaste <- paste('I(', KcatExpression, ' * ', sep = "")
    
    rxnEquations[["elasticity_calc"]] <- rxnEquations[["occupancyEq_list"]]
    rxnEquations[["elasticity_calc"]] <- gsub('(I\\()', Kcatpaste, rxnEquations[["elasticity_calc"]])
    rxnEquations[["elasticity_calc"]] <- sub('I\\(', '\\(', rxnEquations[["elasticity_calc"]])
    
  }
  
  # find the partial derivatives of the kinetic form for each reaction specie
  eq <- eval(parse(text = paste('expression(',rxnEquations[["full_kinetic_form_list"]],')')))
  
  D_full_kinetic_form <- list()
  for(spec in c(kineticPars$rel_spec[(kineticPars$measured & !is.na(kineticPars$measured)) | kineticPars$modelName == "t_metX" & !is.na(kineticPars$modelName)], all_species$rel_spec[all_species$SpeciesType == "Enzyme"])){
    D_full_kinetic_form[[spec]] <- D(eq, spec)
  }
  rxnEquations[["kinetic_form_partials"]] <- D_full_kinetic_form
  
  rxnEquations <<- rxnEquations
  
}  
  
build_kinetic_parameter_priors <- function(rxnSummary, kineticPars, omic_data){
  
  kineticParPrior <- data.frame(distribution = rep(NA, times = nrow(kineticPars)), par_1 = NA, par_2 = NA, par_3 = NA) 
  
  # Options for these parameters are:
  # -@-@ unif: uniform in log-space: par_1 = lowerbound, par_2 = upperbound. draw in log2 space and exponentiate back to linear space
  # -@-@ norm: lognormal: in log2 space draw a value from mean = par_1, sd = par_2
  # -@-@ SpSl: spike and slab (In log2 space: the spike is a point mass at zero with p = par_3, the slab is a normal with mean = par_1 and sd = par_2)
  
  # Specify prior for michaelis constants
  kineticParPrior$distribution[kineticPars$SpeciesType %in% c("Metabolite", "keq")] <- "unif"
  kineticParPrior$par_1[kineticParPrior$distribution == "unif"] <- -15; kineticParPrior$par_2[kineticParPrior$distribution == "unif"] <- 15 # default value to 2^-15:2^15
  
  for(exp_param in kineticPars$formulaName[!is.na(kineticPars$measured) & kineticPars$measured == TRUE]){
    kineticParPrior[kineticPars$formulaName == exp_param & !is.na(kineticPars$modelName), c(2:3)] <- median(omic_data$met_abund[,colnames(omic_data$met_abund) == kineticPars$modelName[kineticPars$formulaName == exp_param & !is.na(kineticPars$formulaName)]]) + c(-15,15)
  }#priors for measured metabolites (either in absolute or relative space) are drawn about the median
  
  # Prior for keq is centered around sum log(sub) - sum log(prod) - this deals with some species being in absolute space and others being absolute measurements
  n_c <- nrow(omic_data$flux)
  kineticParPrior[kineticPars$SpeciesType == "keq", c(2:3)] <- median(rowSums(omic_data$met_abund * c(rep(1,n_c)) %*% t(rxnSummary$rxnStoi[data.table::chmatch(colnames(omic_data$met_abund), names(rxnSummary$rxnStoi))]))) + c(-20,20)
  
  # Specify prior for hill constants
  kineticParPrior$distribution[kineticPars$SpeciesType == "hillCoefficient"] <- "unif"
  kineticParPrior$par_1[kineticPars$SpeciesType == "hillCoefficient"] <- -3
  kineticParPrior$par_2[kineticPars$SpeciesType == "hillCoefficient"] <- 3
  
  return(kineticParPrior)
  
}


#### Fitting reaction equations

fit_reaction_equation_MCMC_NNLS <- function(markov_pars, kineticPars, kineticParPrior, rxnEquations, all_species, omic_data, kinetically_differing_isoenzymes){
  
  ### Initialize parameters & setup tracking of likelihood and parameters ###
  
  lik_track <- NULL
  markov_par_vals <- matrix(NA, ncol = nrow(kineticParPrior), nrow = markov_pars$n_samples)
  colnames(markov_par_vals) <- kineticPars$formulaName
  
  current_pars <- rep(NA, times = nrow(kineticParPrior))
  current_pars <- par_draw(1:nrow(kineticParPrior), current_pars, kineticParPrior)
  
  current_lik <- lik_calc_fittedSD(current_pars, kineticPars, all_species, rxnEquations, omic_data, kinetically_differing_isoenzymes)
  
  ### Generate markov chain ###
  
  for(i in 1:(markov_pars$burn_in + markov_pars$n_samples*markov_pars$sample_freq)){
    for(j in 1:nrow(kineticPars)){#loop over parameters values
      proposed_par <- par_draw(j, current_pars, kineticParPrior)
      proposed_lik <- lik_calc_fittedSD(proposed_par, kineticPars, all_species, rxnEquations, omic_data, kinetically_differing_isoenzymes)
      if(runif(1, 0, 1) < exp(proposed_lik - current_lik) | (proposed_lik == current_lik & proposed_lik == -Inf)){
        current_pars <- proposed_par
        current_lik <- proposed_lik
      }
    }
    
    if(i > markov_pars$burn_in){
      if((i - markov_pars$burn_in) %% markov_pars$sample_freq == 0){
        markov_par_vals[(i - markov_pars$burn_in)/markov_pars$sample_freq,] <- current_pars
        lik_track <- c(lik_track, current_lik)
      }
    }
  }
  
  colnames(markov_par_vals) <- kineticPars$rel_spec
  
  markov_par_vals <<- markov_par_vals  
  lik_track <<- lik_track
  
  }

par_draw <- function(updates, current_pars, kineticParPrior){
  ### Update parameters using their prior (given by kineticParPrior) - update those those parameters whose index is in "updates" ###
  ### Parameters are all returned in log-space (base e) ###
  
  draw <- current_pars
  for(par_n in updates){
    if(kineticParPrior$distribution[par_n] == "unif"){
      draw[par_n] <- runif(1, kineticParPrior$par_1[par_n], kineticParPrior$par_2[par_n])
    } else if(kineticParPrior$distribution[par_n] == "norm"){
      draw[par_n] <- rnorm(1, kineticParPrior$par_1[par_n], kineticParPrior$par_2[par_n])
    } else if(kineticParPrior$distribution[par_n] == "SpSl"){
      draw[par_n] <- ifelse(rbinom(1, 1, kineticParPrior$par_3[par_n]) == 0, 0, rnorm(1, kineticParPrior$par_1[par_n], kineticParPrior$par_2[par_n]))
    } else {
      stop("invalid distribution")
    }
  }
  draw
}
 
lik_calc_fittedSD <- function(proposed_params, kineticPars, all_species, rxnEquations, omic_data, kinetically_differing_isoenzymes){

  #### determine the likelihood of predicted flux as a function of metabolite abundance and kinetics parameters relative to actual flux ####
  
  n_c <- nrow(omic_data$met_abund)
  par_stack <- rep(1, n_c) %*% t(proposed_params); colnames(par_stack) <- kineticPars$formulaName
  
  par_stack <- 2^par_stack
  occupancy_vals <- data.frame(omic_data$met_abund, par_stack)
  
  if(!(kinetically_differing_isoenzymes)){
    predOcc <- eval(rxnEquations[["l_occupancyExpression"]], occupancy_vals) #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
    enzyme_activity <- (predOcc %*% t(rep(1, sum(all_species$SpeciesType == "Enzyme"))))*2^omic_data$enzyme_abund #occupany of enzymes * relative abundance of enzymes
  }else{
    enzyme_activity <- NULL
    for(isoenzyme in names(rxnSummary$rxnForm)){
      predOcc <- eval(rxnEquations[["l_occupancyExpression"]][[isoenzyme]], occupancy_vals)
      enzyme_activity <- cbind(enzyme_activity, predOcc %*% t(rep(1, sum(occEqtn_complex_match$occEqtn == isoenzyme)))*2^omic_data$enzyme_abund[,colnames(omic_data$enzyme_abund) %in% occEqtn_complex_match$complex[occEqtn_complex_match$occEqtn == isoenzyme]])
    }
  }
  
  # fit flux ~ enzyme*occupancy using non-negative least squares (all enzymes have activity > 0, though negative flux can occur through occupancy)
  # flux objective is set as the average of the minimal and maximal allowable flux flowing through the reaction at the optimal solution
  
  flux_fit <- nnls(enzyme_activity, (omic_data$flux$FVAmax + omic_data$flux$FVAmin)/2) 
  
  # setting the variance from residual mean-squared error
  
  nparam <- sum(kineticPars$measured[kineticPars$SpeciesType == "Metabolite"]) + # number of measured metabolites
    sum(all_species$SpeciesType == "Enzyme") + # number of measured enzyme groups
    sum(kineticPars$SpeciesType == "keq") # plus 1 if keq is included
  fit_resid_error <- sqrt(mean((flux_fit$resid - mean(flux_fit$resid))^2))*sqrt(n_c/(n_c-nparam))
  
  # integrate over the cdf from FVA_min to FVA_max t
  lik <- (pnorm(omic_data$flux$FVAmax, flux_fit$fitted, fit_resid_error) - pnorm(omic_data$flux$FVAmin, flux_fit$fitted, fit_resid_error))/(omic_data$flux$FVAmax - omic_data$flux$FVAmin)
  
  return(sum(log(lik)))
  
  }

#### Misc

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

