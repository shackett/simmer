options(stringsAsFactors = F)

###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@
##### Load files describing valid reactions, species (their composition) both from the core SBML model and supplemented manual annotations #####
###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@

load_metabolic_model()

###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@##
##### Load files describing boundary conditions, reaction reversibility and auxotrophies #####
###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@##

format_boundary_conditions()

###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@
####### Setup matrices defining the stoichiometry of each reaction and how reactions will be constrained #######
###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@

# identify reaction that cannot carry flux based on inputs, outputs and directionality
infRxMet <- setup_FBA_constraints(return_infeasible = T, infRxMet = NULL)

setup_FBA_constraints(return_infeasible = F, infRxMet)

inferred_fluxes <- calculate_QP_fluxes()

read_write_python_FVA(status = "write", directory = "~/Desktop/")



#### Setting up stoichiometric inputs

load_metabolic_model <- function(){
  
  ### Load core model derived from SBML reconstruction
  
  # reaction stoichiometry, compartment, metabolite / rxn IDs
  rxnFile = read.delim("companionFiles/genome_scale_model/rxn_yeast.tsv")
  # reaction enzymes and annotations
  rxnparFile = read.delim("companionFiles/genome_scale_model/rxn_par_yeast.tsv")
  # mapping between unique metabolites and compartmentalized metabolites
  corrFile = read.delim("companionFiles/genome_scale_model/spec_yeast.tsv")
  # metabolite annotations (ChEBI and KEGG)
  specparFile = read.delim("companionFiles/genome_scale_model/species_par_yeast.tsv")
  # comparment IDs to compartment names
  compFile = read.delim("companionFiles/genome_scale_model/comp_yeast.tsv")
  # model-specified reaction directionality
  fluxDirFile = read.delim("companionFiles/genome_scale_model/flux_dir_yeast.tsv")
  
  ### Add additional reactions ###
  
  customList <- parse_custom("companionFiles/genome_scale_model/customRxns.txt")
  
  rxnFile <- rbind(rxnFile, customList$rxnFile)
  rxnparFile <- rbind(rxnparFile, customList$rxnparFile)
  corrFile <- rbind(corrFile, customList$corrFile)
  specparFile <- rbind(specparFile, customList$specparFile)
  fluxDirFile <- rbind(fluxDirFile, customList$fluxDirFile)
  
  assign("rxnFile", rxnFile, envir=globalenv())
  assign("rxnparFile", rxnparFile, envir=globalenv())
  assign("corrFile", corrFile, envir=globalenv())
  assign("specparFile", specparFile, envir=globalenv())
  assign("fluxDirFile", fluxDirFile, envir=globalenv())
  
  ### Determine unique metabolites and reactions ###
  
  reactions <- sort(unique(rxnFile$ReactionID))
  rxnStoi <- rxnFile[is.na(rxnFile$StoiCoef) == FALSE,]
  metabolites <- sort(unique(rxnStoi$Metabolite))
  
  assign("reactions", reactions, envir=globalenv())
  assign("metabolites", metabolites, envir=globalenv())
  
  # generate the stoichiometric matrix from reactions
  
  stoiMat <- build_stoiMat(metabolites, reactions, corrFile, rxnFile, internal_names = TRUE)
  assign("stoiMat", stoiMat, envir=globalenv())
  
  # generated a named version of the stoichiometric matrix (this is used to search for some reactions)
  
  named_stoi <- stoiMat
  met_dict <- metIDtoSpec(rownames(named_stoi)); met_dict <- sapply(c(1:length(named_stoi[,1])), function(x){met_dict[x][[1]]})
  rxn_dict <- rxnIDtoEnz(colnames(named_stoi)); rxn_dict <- sapply(c(1:length(named_stoi[1,])), function(x){rxn_dict[x][[1]]})
  rownames(named_stoi) <- met_dict; colnames(named_stoi) <- rxn_dict
  assign("named_stoi", named_stoi, envir=globalenv())
  
  # Determine the compartmentation of each reaction

  compartment <- sapply(reactions, function(x){rxnFile$Compartment[rxnFile$ReactionID == x][1]})
  assign("compartment", compartment, envir=globalenv())
  
  ### update reaction directionalities
  
  reversibleRx <- data.frame(rx = reactions, reversible = 0, CCdG = NA, CCdGsd = NA, CCdGdir = NA, modelRev = NA, modelBound = NA, manual = NA, rxFlip = NA, annotComment = NA)

  reversibleRx$modelRev = fluxDirFile$Reversible[data.table::chmatch(reversibleRx$rx, fluxDirFile$ReactionID)]
  reversibleRx$modelBound = fluxDirFile$FluxBound[data.table::chmatch(reversibleRx$rx, fluxDirFile$ReactionID)]
  reversibleRx$reversible = ifelse(reversibleRx$modelBound == "greaterEqual", 1, 0) # directionality is coded as -1: irreversibly backward, 0: reversible, 1: irreversibly forward

  # append directionality with manual annotation in several cases
  
  manualDirectionality <- read.delim("companionFiles/genome_scale_model/thermoAnnotate.txt")
  reversibleRx$manual[data.table::chmatch(manualDirectionality$Reaction, reversibleRx$rx)] <- manualDirectionality$Direction
  reversibleRx$reversible[!is.na(reversibleRx$manual)] <- reversibleRx$manual[!is.na(reversibleRx$manual)]
  
  assign("reversibleRx", reversibleRx, envir=globalenv())
  
  ### Associate reactions with pathway and measured proteins to favor measured reactions and major pathways
  
  # Associate proteins with reactions and format complexes
  
  rxn_enzyme_groups = NULL
  for(rxN in 1:nrow(rxnparFile)){
    
    enzymeCombos <- rxnparFile$Enzymes[rxN]
    
    if(enzymeCombos == ""){
      next
    } # no enzymes
    
    enzymeCombos <- strsplit(enzymeCombos, ' OR ')[[1]] # split complexes demarcated by OR
    
    for(combo in 1:length(enzymeCombos)){
      rxn_enzyme_groups <- rbind(rxn_enzyme_groups,
                                 data.frame(reaction = rxnparFile$ReactionID[rxN],
                                            group = combo,
                                            enzyme = regmatches(enzymeCombos[combo], gregexpr('[YQ][A-Z0-9]+[WC]?-?[A-Z]{0,1}', enzymeCombos[combo]))[[1]]))
    }
  }

  enzyme_abund <- read.delim("companionFiles/flux_input_data/proteinAbundance.tsv")
  rownames(enzyme_abund) <- enzyme_abund$Gene; enzyme_abund <- enzyme_abund[,-1]

  prot_matches <- sapply(reactions, function(x){
    rxMatches <- rxn_enzyme_groups$enzyme[rxn_enzyme_groups$reaction == x]
    length(rownames(enzyme_abund)[rownames(enzyme_abund) %in% rxMatches]) != 0
  })

  # Determine which metabolic pathways are associated with a reaction ###

  reactionMatches <- data.frame(reactionID = rxnparFile$ReactionID[grep('reaction/', rxnparFile$Annotation)],
                                KEGGrxnID = sapply(grep('reaction/', rxnparFile$Annotation, value = T), function(x){regmatches(x, regexpr('R[0-9]+', x))}))

  reactions_to_pathways = read.delim("http://rest.kegg.jp/link/reaction/pathway", header = FALSE); colnames(reactions_to_pathways) <- c("pathwayCode", "RID")
  reactions_to_pathways$RID = sapply(reactions_to_pathways$RID, function(x){sub('rn:', '', x)})

  rxPathways = read.delim("http://rest.kegg.jp/list/pathway", header = FALSE); colnames(rxPathways) = c("pathwayCode", "pathway")
  reactions_to_pathways$pathway = rxPathways$pathway[data.table::chmatch(reactions_to_pathways$pathwayCode, rxPathways$pathwayCode)]

  reactionMatches$pathway = sapply(reactionMatches$KEGGrxnID, function(x){
    PS = reactions_to_pathways$pathway[reactions_to_pathways$RID == x]
    PS <- PS[!is.na(PS)]
    paste(PS, collapse = "__")
  })
  
  # favor flux through glycolysis and TCA

  pathways <- sapply(reactionMatches$pathway, function(x){strsplit(x, '__')[[1]]})
  unq_pathways <- unique(unlist(pathways))
  centralC <- c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)")

  centralCmatch <- unlist(lapply(pathways, function(x){
    sum(x %in% centralC) != 0
  }))

  centralCmatch[reactionMatches$reactionID %in% c("r_0713", "r_0962")] <- FALSE # remove MDH and PyK from this list to remove futile cycling

  centralCrxnMatch <- rep(FALSE, times = length(reactions))
  centralCrxnMatch[reactions %in% reactionMatches$reactionID[centralCmatch]] <- TRUE

  ETCrxns <- c("r_0226", "r_0438", "r_0439", "r_0773", "r_1021", "r_0831", "r_0832") # electron transport reactions

  centralCrxnMatch[reactions %in% ETCrxns] <- TRUE

  # favor flux through cytosolic ATPase, + transport of ATP, ADP + Pi to restrict the wasting of excess energy to a single reaction

  ATPbreakdownRxns <- c("r_4042", "r_1110", "r_1111", "r_1661", "r_3543", "r_3585", "r_3601", "r_3651", "r_3666",
                      "r_1244", "r_1245", "r_2005", "r_2008", "r_3537", "r_3605", "r_3649", "r_3663", "r_3940", "r_3961")

  prot_penalty <- (1 - prot_matches)/2 + (1 - centralCrxnMatch)/2 # penalization by fraction of non-measured enzymes and favor central C metabolism
  prot_penalty[reactions %in% ATPbreakdownRxns] <- 0.1

  freeTransportRxns = c("r_1277", "r_2096", "r_1978", "r_1979", "r_1696", "r_1697") # transport of water, carbon dioxide and oxygen
  prot_penalty[reactions %in% freeTransportRxns] <- 0
  
  assign("prot_penalty", prot_penalty, envir=globalenv())
  
}

parse_custom <- function(customRx){
  
  ### add additional reactions, metabolites, etc to metabolic model
  
  outputList <- list()
  
  inputFile = read.table(customRx, header = F, sep = "\t", fill = T, blank.lines.skip = F)
  
  ### Species-level annotation ###
  
  input_bound = c(1:nrow(inputFile))[inputFile[,1] == "!Species"] + c(1,-1)
  spec_input = inputFile[(input_bound[1]+1):input_bound[2],colSums(inputFile[(input_bound[1]+1):input_bound[2],] != "") != 0]
  colnames(spec_input) <- inputFile[input_bound[1],colSums(inputFile[(input_bound[1]+1):input_bound[2],] != "") != 0]
  
  corr = spec_input[,colnames(spec_input) %in% c("SpeciesID", "SpeciesName", "SpeciesType", "Compartment")]
  outputList$corrFile = corr[match(unique(corr$SpeciesType), corr$SpeciesType),]

  spec_par = spec_input[,colnames(spec_input) %in% c("SpeciesID", "SpeciesName", "Annotation")]
  outputList$specparFile = spec_par
  
  ### Reaction stoichiometry annotation ###
  
  input_bound = c(1:nrow(inputFile))[inputFile[,1] == "!Reactions"] + c(1,-1)
  rxn_input = inputFile[(input_bound[1]+1):input_bound[2],]
  colnames(rxn_input) <- inputFile[input_bound[1],]
  
  outputList$rxnFile = rxn_input
  
  ### Reaction flux and references annotation ###
  
  input_bound = c(1:nrow(inputFile))[inputFile[,1] == "!Reaction_Parameters"] + c(1,-1)
  rxPar_input = inputFile[(input_bound[1]+1):input_bound[2],apply(inputFile[(input_bound[1]+1):input_bound[2],], 2, function(x){sum(x[!is.na(x)] != "")}) != 0]
  colnames(rxPar_input) <- inputFile[input_bound[1],apply(inputFile[(input_bound[1]+1):input_bound[2],], 2, function(x){sum(x[!is.na(x)] != "")}) != 0]
  
  outputList$rxnparFile = rxPar_input[,colnames(rxPar_input) %in% c("ReactionID", "Enzymes", "Annotation")]
  outputList$fluxDirFile = rxPar_input[,colnames(rxPar_input) %in% c("ReactionID", "Reversible", "FluxBound")]
  
  return(outputList)
  
}

build_stoiMat <- function(metabolites, reactions, corrFile, rxnFile, internal_names = FALSE){
  
  # create the stoichiometry matrix, indicating the change of metabolites (rows) per chemical reaction (columns)
  
  stoiMat <- matrix(data = 0, ncol = length(reactions), nrow = length(metabolites))
  
  metName = rep(NA, times = length(metabolites))
  for (i in 1:length(metabolites)){
    metName[i] = as.character(corrFile$SpeciesName[corrFile$SpeciesID %in% metabolites[i]][1])
  }
  rxnName = rep(NA, times = length(reactions))
  for (i in 1:length(reactions)){
    rxnName[i] = as.character(rxnFile$Reaction[rxnFile$ReactionID %in% reactions[i]][1])
  }
  rownames(stoiMat) <- metName	
  colnames(stoiMat) <- rxnName
  
  rxStoi <- rxnFile[!is.na(rxnFile$StoiCoef),]
  for (i in 1:nrow(rxStoi)){
    stoiMat[c(1:length(metabolites))[metabolites == rxStoi$Metabolite[i]], c(1:length(reactions))[reactions == rxStoi$ReactionID[i]]] <- as.numeric(rxStoi$StoiCoef[i])
  }
  
  if (internal_names == TRUE){
    rownames(stoiMat) <- metabolites	
    colnames(stoiMat) <- reactions
  }
  
  return(stoiMat)
  
}


#### Setting up boundary fluxes

format_boundary_conditions <- function(){
  
  require(dplyr)
  
  ### Format experimental flux inputs and outputs for each condition (treatment_par) and match species to model reconstruction
  # load experimental data
  # combine boundary fluxes and uncertainties for each treatment into a list: treatment_par
  # flag the entries in the metabolic model that match boundary species
  
  ### Experimental inputs (media formulation, uptake, excretion, incorporation rates)
  
  nutrientFile <- read.delim("companionFiles/flux_input_data/Boer_nutrients.txt")[1:6,1:6]
  nutrientCode <- data.frame(nutrient = colnames(nutrientFile)[-1], shorthand = c("N", "P", "C", "L", "U"))
  
  rownames(nutrientFile) <- nutrientFile[,1]; nutrientFile <- nutrientFile[,-1]
  nutModelNames <- data.frame(commonName = rownames(nutrientFile), modelName = sapply(rownames(nutrientFile), function(x){paste(x, '[extracellular]')}))
  rownames(nutrientFile) <- nutModelNames$modelName

  load("companionFiles/flux_input_data/boundaryFluxes.Rdata") #load condition specific boundary fluxes and chemostat info (actual culture DR)

  
  ### Define the treatment in terms of nutrient availability and auxotrophies

  treatment_par <- list()
  n_c <- nrow(chemostatInfo)
  for(i in 1:n_c){
    #define nutrient uptake and excretion rate - soft matches on (using maximal available for now)
    measured_bounds <- mediaSummary[mediaSummary$condition == chemostatInfo$ChemostatCond[i],]
    measured_bounds <- rbind(measured_bounds, data.frame(condition = chemostatInfo$ChemostatCond[i], specie = rownames(nutrientFile)[!(rownames(nutrientFile) %in% measured_bounds$specie)], change = NA, sd = NA, lb = 0, 
                                                         ub = nutrientFile[,colnames(nutrientFile) == nutrientCode$nutrient[nutrientCode$shorthand == chemostatInfo$Limitation[i]]][!(rownames(nutrientFile) %in% measured_bounds$specie)], type = "uptake", density = NA))
    
    # multiply steady-state concentrations by DR to get the uptake/excretion rates
    measured_bounds$change <- measured_bounds$change*chemostatInfo$actualDR[i]
    measured_bounds$sd <- measured_bounds$sd*chemostatInfo$actualDR[i]
    measured_bounds$lb <- measured_bounds$lb*chemostatInfo$actualDR[i]
    measured_bounds$ub <- measured_bounds$ub*chemostatInfo$actualDR[i]
    
    # remove phosphate because empirical uptake rates far exceed capacity of biomass assimilation
    measured_bounds <- data.frame(measured_bounds)
    measured_bounds[measured_bounds$specie == "phosphate [extracellular]", colnames(measured_bounds) %in% c("change", "sd")] <- NA
    measured_bounds <- data.table::data.table(measured_bounds)
    
    # define approximate oxygen uptake using a RQ of 5 for non-carbon-limited culture and < 5 for carbon-limited culture
    # employed by approximating vCO2 as 5/4 [ethanol + actetate]
    
    if(chemostatInfo$Limitation[i] == "C"){
      oxygen_uptake = (measured_bounds$change[measured_bounds$specie %in% c("D-glucose [extracellular]")]/5 + 
                         sum(measured_bounds$change[measured_bounds$specie %in% c("ethanol [extracellular]", "acetate [extracellular]")])/2)/2    
      
      oxygen_bounds <- data.frame(condition = chemostatInfo$ChemostatCond[i], specie = 'oxygen [extracellular]', change = 0, sd = Inf, 
                                  lb = oxygen_uptake,
                                  ub = Inf, type = "uptake", density = NA)
    } else {
      oxygen_uptake = (measured_bounds$change[measured_bounds$specie %in% c("D-glucose [extracellular]")]/5 + 
                         sum(measured_bounds$change[measured_bounds$specie %in% c("ethanol [extracellular]", "acetate [extracellular]")])/4)/2    
      
      oxygen_bounds <- data.frame(condition = chemostatInfo$ChemostatCond[i], specie = 'oxygen [extracellular]', 
                                  change = oxygen_uptake,
                                  sd = NA, lb = 0, ub = Inf, type = "uptake", density = NA)
      oxygen_bounds$sd <- oxygen_bounds$change/5
    }
    
    measured_bounds <- rbind(measured_bounds, oxygen_bounds)
    
    treatment_par[[chemostatInfo$ChemostatCond[i]]][["nutrients"]] <- measured_bounds
    
    #define ura3 and leu2 auxotrophies
    if(chemostatInfo$Limitation[i] == "L"){treatment_par[[chemostatInfo$ChemostatCond[i]]][["auxotrophies"]] <- as.character(unique(rxnFile[grep("isopropylmalate dehydrogenase", rxnFile$Reaction),]$ReactionID))}
    if(chemostatInfo$Limitation[i] == "U"){treatment_par[[chemostatInfo$ChemostatCond[i]]][["auxotrophies"]] <- as.character(unique(rxnFile[grep("orotidine", rxnFile$Reaction),]$ReactionID))}
    if(chemostatInfo$Limitation[i] %in% c("C", "P", "N")){treatment_par[[chemostatInfo$ChemostatCond[i]]][["auxotrophies"]] <- NA}
    
    #define observed fluxes per culture volume #scaled to intracellular volume later
    biomass_match <- data.frame(specie = comp_by_cond$compositionFile$MetName, AltName = comp_by_cond$compositionFile$AltName,change = unname(-1*comp_by_cond$cultureMolarity[,colnames(comp_by_cond$cultureMolarity) == chemostatInfo$ChemostatCond[i]]*chemostatInfo$actualDR[i]), index = comp_by_cond$compositionFile$index)
    biomass_list <- list()
    
    for(component in unique(comp_by_cond$compositionFile$varCategory)){
      principal_costs <- biomass_match[comp_by_cond$compositionFile$varCategory %in% component,]
      
      if(component == "Maintenance ATP hydrolysis"){
        total_costs <- principal_costs[,colnames(principal_costs) != "index"]
      }else{
        #costs of monomer assimilation incorporated into biomass flux
        energetic_costs <- as.matrix(comp_by_cond$biomassExtensionE[principal_costs$index,-1])
        energetic_costs_aggregate <- t(principal_costs$change) %*% energetic_costs; colnames(energetic_costs_aggregate) <- colnames(comp_by_cond$biomassExtensionE)[-1]
        total_costs <- rbind(principal_costs[,colnames(principal_costs) != "index"], data.frame(specie = colnames(energetic_costs_aggregate), AltName = colnames(energetic_costs_aggregate), change = t(unname(energetic_costs_aggregate)))[energetic_costs_aggregate != 0,])
      }
      biomass_list[[component]]$exchange = total_costs
      
      # define the accuracy of a constraint in terms of the coefficient of variation - sd over mean
      if(component %in% colnames(comp_by_cond$CV_table)){
        biomass_list[[component]]$SD = as.numeric(subset(comp_by_cond$CV_table, comp_by_cond$CV_table$condition == chemostatInfo$ChemostatCond[i], component))
      }else{
        biomass_list[[component]]$SD = 1/10
      }
    }
    
    treatment_par[[chemostatInfo$ChemostatCond[i]]][["boundaryFlux"]] = biomass_list
  }
  
  assign("treatment_par", treatment_par, envir=globalenv())
  assign("comp_by_cond", comp_by_cond, envir=globalenv())
  assign("chemostatInfo", chemostatInfo, envir=globalenv())
  
  possibleAuxotrophies <- c(as.character(unique(rxnFile[grep("isopropylmalate dehydrogenase", rxnFile$Reaction),]$ReactionID)), as.character(unique(rxnFile[grep("orotidine", rxnFile$Reaction),]$ReactionID)))
  assign("possibleAuxotrophies", possibleAuxotrophies, envir=globalenv())
  
  ### Define species involved in boundary-conditions
  
  # extract the metabolite ID corresponding to the extracellular introduction of nutrients
  boundary_met <- treatment_par[[1]]$nutrients %>% filter(type == "uptake") %>% 
    select(SpeciesName = specie) %>% tbl_df() %>% left_join(corrFile, by = "SpeciesName")
  if(nrow(boundary_met) != nrow(treatment_par[[1]]$nutrients %>% filter(type == "uptake"))){stop("A valid match was not found for all absorbed metabolites")}
  assign("boundary_met", boundary_met, envir=globalenv())

  # extract the IDs of excreted metabolites
  excreted_met <- treatment_par[[1]]$nutrients %>% filter(type == "excretion") %>% 
    select(SpeciesName = specie) %>% tbl_df() %>% left_join(corrFile, by = "SpeciesName")
  if(nrow(excreted_met) != nrow(treatment_par[[1]]$nutrients %>% filter(type == "excretion"))){stop("A valid match was not found for all excreted metabolites")}
  assign("excreted_met", excreted_met, envir=globalenv())
  
  # extract the metabolite ID corresponding to cytosolic metabolites being assimilated into biomass
  sinks <- unique(c(comp_by_cond$compositionFile$AltName[comp_by_cond$compositionFile$FluxType == "Boundary"], colnames(comp_by_cond$biomassExtensionE)[-1]))
  comp_met <- data.frame(SpeciesName = sinks) %>% tbl_df() %>% left_join(corrFile, by = "SpeciesName")
  if(nrow(comp_met) != length(sinks)){stop("A valid match was not found for all sinks")}
  assign("comp_met", comp_met, envir=globalenv())
  
  # freely exchanging metabolites through extracellular compartment ##
  free_flux <- c("carbon dioxide [extracellular]", "H2O [extracellular]", "H+ [extracellular]")
  freeExchange_met <- data.frame(SpeciesName = free_flux) %>% tbl_df() %>% left_join(corrFile, by = "SpeciesName")
  assign("freeExchange_met", freeExchange_met, envir=globalenv())
  
}


#### Format optimization problem

setup_FBA_constraints <- function(return_infeasible = F, infRxMet = NULL){
  
  # if return_infeasible is TRUE then a data.frame of reactions that cannot carry flux (and associated metabolites) will be returned
  # providing this data.frame as infRxMet will exclude these reactions and species from subsequent steps
  
  stopifnot(is.logical(return_infeasible),
            class(infRxMet) %in% c("NULL", "character"))
  
  ### Remove reactions which are incompatable with our method - e.g. an invariant biomass function
  
  labelz <- c("biomass")
  aggregate_rxns <- NULL
  for(l in length(labelz)){
    aggregate_rxns <- union(aggregate_rxns, rxn_search(labelz[l], named_stoi, is_rxn = TRUE, index = TRUE))
  }

  # Remove lipid breakdown reactions
  # 1) reactions are annotated as a hydrolase or lipase.
  # 2) reactions produce a fatty acid annotated with the isa fatty acid tag

  fatty_acids <- data.frame(name = c("myristate", "palmitate", "palmitoleate", "stearate", "oleate"), tID = NA)
  fatty_acids$tID <- sapply(fatty_acids$name, function(x){unique(corrFile$SpeciesType[grep(paste('^', x, sep = ""), corrFile$SpeciesName)])})
  all_FA <- corrFile$SpeciesID[corrFile$SpeciesType %in% fatty_acids$tID]

  HLsubset <- stoiMat[metabolites %in% all_FA, grep('hydrolase|lipase', colnames(named_stoi))]
  lipase_reactions <- colnames(HLsubset)[colSums(HLsubset) != 0]
  
  # remove the above reactions from the model
  
  S_rxns <- stoiMat[,!(colnames(stoiMat) %in% c(aggregate_rxns, lipase_reactions))]
  
  if(!is.null(infRxMet)){
    # if infRxMet are included, remove these reactions as well
    infRxMet <- data.frame(infRxMet); colnames(infRxMet) <- "ID"
    infRxMet <- rbind(infRxMet, 'r_0163') # remove several reactions which wont be used and cycle (ADH reverse)
  
    infRxMet$type <- substr(infRxMet$ID, 1, 1)
  
    reversibleRx <- reversibleRx[!(reversibleRx$rx %in% infRxMet$ID[infRxMet$type == "r"]),]
    metabolites <- metabolites[!(rownames(S_rxns) %in% infRxMet$ID[infRxMet$type == "s"])]
    
    assign("reversibleRx", reversibleRx, envir=globalenv())
    assign("metabolites", metabolites, envir=globalenv())
    
    S_rxns <- S_rxns[!(rownames(S_rxns) %in% infRxMet$ID[infRxMet$type == "s"]), !(colnames(S_rxns) %in% infRxMet$ID[infRxMet$type == "r"])]
  }
  
  assign("S_rxns", S_rxns, envir=globalenv())
  
  ### Set up the equality and inequality constriants for FBA
  
  setup_equality_constraints()
  
  setup_inequality_constraints()
  
  ### F : flux balance

  Fzero <- rep(0, times = length(S[,1]))
  assign("Fzero", Fzero, envir=globalenv())
  
  ### Identify reactions that cannot carry flux based on inputs and outputs
  
  if(return_infeasible == TRUE){
    
    S_tmp <- S # S_tmp needs to have temporary consumption rates filled in
    S_tmp_blank <- S_tmp[,colnames(S_tmp) %in% names(biomassConv)]
    for(j in 1:ncol(S_tmp_blank)){
      biomass_fill <- biomassConv[[colnames(S_tmp_blank)[j]]]$conversion
      biomass_fill$coef <- treatment_par[[1]]$boundaryFlux[[strsplit(colnames(S_tmp_blank)[j], split = ' comp')[[1]][1]]]$exchange$change
      S_tmp_blank[biomass_fill$index,j] <- biomass_fill$coef  
    }
    
    S_tmp_blank <- S_tmp_blank *  rep(1, nrow(S_tmp)) %*% t(ifelse(Sinfo$direction[colnames(S_tmp) %in% names(biomassConv)] == "F", 1, -1))
    
    S_tmp[,colnames(S_tmp) %in% names(biomassConv)] <- S_tmp_blank
    
    cpdfil_old = rep(T,nrow(S_tmp))
    cpdfil = rep(F,nrow(S_tmp))
    while (any(cpdfil != cpdfil_old)){
      cpdfil_old = cpdfil
      cpdfil = apply(S_tmp,1,function(x){all(x>=0)||all(x<=0)}) #metabolites which are only consumed or produced
      rxnfil = apply(S_tmp[cpdfil,],2,function(x){any(x!=0)}) #reactions involving these metabolites are flagged
      S_tmp[,rxnfil]=0 #flagged reactions are zeroed
    }
    rxnfil = apply(S_tmp,2,function(x){all(x==0)})
    infRxnMet = c(rownames(S_tmp)[cpdfil],unique(substr(colnames(S_tmp)[rxnfil & grepl('^r_',colnames(S_tmp))],1,6)))
    
    return(infRxnMet)
    
  }
}

setup_equality_constraints <- function(){
  
  ### Split all reacitons into a positive and a negative component
  
  validrxns <- reversibleRx[reversibleRx$rx %in% colnames(S_rxns),]
  validrxnsFluxdir <- sapply(colnames(S_rxns), function(rxdir){
    validrxns$reversible[validrxns$rx == rxdir]
  })
  
  stoiRxSplit <- NULL
  for(rxsplit in 1:length(validrxnsFluxdir)){
    if(unname(validrxnsFluxdir[rxsplit]) == 0){out <- rbind(c(paste(c(names(validrxnsFluxdir)[rxsplit], "F"), collapse = '_'), names(validrxnsFluxdir)[rxsplit], "F"), c(paste(c(names(validrxnsFluxdir)[rxsplit], "R"), collapse = '_'), names(validrxnsFluxdir)[rxsplit], "R"))}
    if(unname(validrxnsFluxdir[rxsplit]) == 1){out <- c(paste(c(names(validrxnsFluxdir)[rxsplit], "F"), collapse = '_'), names(validrxnsFluxdir)[rxsplit], "F")}  
    if(unname(validrxnsFluxdir[rxsplit]) == -1){out <- c(paste(c(names(validrxnsFluxdir)[rxsplit], "R"), collapse = '_'), names(validrxnsFluxdir)[rxsplit], "R")}  
    stoiRxSplit <- rbind(stoiRxSplit, out)
  }
  stoiRxSplit <- as.data.frame(cbind(stoiRxSplit)); colnames(stoiRxSplit) <- c("rxDesignation", "reaction", "direction")
  
  S_rxns_split <- sapply(stoiRxSplit$reaction, function(rx){
    S_rxns[,colnames(S_rxns) == rx]
  })
  colnames(S_rxns_split) <- stoiRxSplit$rxDesignation
  
  ### Stoichiometry and bounds of boundary fluxes
  
  # Unconstrained Chemical Influx - e.g. water, gases
  
  free_rxns <- sapply(freeExchange_met$SpeciesID, function(id){c(1:length(metabolites))[metabolites == id]})
  
  freeS <- matrix(0, ncol = length(free_rxns), nrow = length(metabolites))
  for(i in 1:length(free_rxns)){
    freeS[free_rxns[i],i] <- 1
  }
  
  freeRxSplit <- data.frame(rxDesignation = c(paste(paste(c(freeExchange_met$SpeciesName), "boundary"), "F", sep = '_'), paste(paste(c(freeExchange_met$SpeciesName), "boundary"), "R", sep = '_')), 
                            reaction = paste(c(freeExchange_met$SpeciesName, freeExchange_met$SpeciesName), "boundary"), direction = rep(c("F", "R"), each = length(freeExchange_met$SpeciesName))
  )
  
  freeS_split <- cbind(freeS, freeS)
  colnames(freeS_split) <- freeRxSplit$rxDesignation
  
  
  # Nutrient Influx
  
  nutrient_rxns <- sapply(boundary_met$SpeciesID, function(id){c(1:length(metabolites))[metabolites == id]})
  nutrientS <- matrix(0, ncol = length(nutrient_rxns), nrow = length(metabolites))
  for(i in 1:length(nutrient_rxns)){
    nutrientS[nutrient_rxns[i],i] <- 1
  }
  
  # oxygen is directly coupled with ETC, to avoid un-physiological cycling ##
  oxygen_stoi <- S_rxns[,grep('r_0438', colnames(S_rxns))]
  oxygen_stoi[names(oxygen_stoi) == "s_1278"] <- 0
  nutrientS[,boundary_met$SpeciesName == "oxygen [extracellular]"] <- oxygen_stoi # oxygen uptake == complex IV flux
  S_rxns_split[grep('s_2817', rownames(S_rxns_split)), grep('r_218[23]', colnames(S_rxns_split))] <- 0 # make oxygen use for FA-desaturation a freebee.
  
  nutrientRxSplit <- data.frame(rxDesignation = c(paste(c(boundary_met$SpeciesName), "boundary_offset"), paste(c(boundary_met$SpeciesName), "boundary_match_F"), paste(c(boundary_met$SpeciesName), "boundary_match_R")), 
                                reaction = rep(paste(c(boundary_met$SpeciesName), "boundary"), times = 3), direction = c(rep("F", length(nutrient_rxns)*2), rep("R", length(nutrient_rxns))))
  
  nutrientS_split <- cbind(nutrientS, nutrientS, nutrientS)
  colnames(nutrientS_split) <- nutrientRxSplit$rxDesignation
  
  
  ## Excreted Metabolite Efflux - non-negative
  
  efflux_rxns <- sapply(excreted_met$SpeciesID, function(id){c(1:length(metabolites))[metabolites == id]})
  effluxS <- matrix(0, ncol = length(efflux_rxns), nrow = length(metabolites))
  for(i in 1:length(efflux_rxns)){
    effluxS[efflux_rxns[i],i] <- -1
  }
  
  effluxRxSplit <- data.frame(rxDesignation = c(paste(c(excreted_met$SpeciesName), "boundary_offset"), paste(c(excreted_met$SpeciesName), "boundary_match_F"), paste(c(excreted_met$SpeciesName), "boundary_match_R")), 
                              reaction = rep(paste((excreted_met$SpeciesName), "boundary"), times = 3), direction = c(rep("F", length(efflux_rxns)*2), rep("R", length(efflux_rxns))))
  
  effluxS_split <- cbind(effluxS, effluxS, effluxS)
  colnames(effluxS_split) <- effluxRxSplit$rxDesignation
  
  
  ## Composition fxn ##
  
  ## generate a conversion matrix between each variance categories components and their indecies ##
  ## the stoichiometery will be populated on a by-condition basis ##
  
  biomassS <- matrix(0, ncol = length(unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Boundary"])), nrow = length(metabolites))
  
  biomassRxSplit <- data.frame(rxDesignation = c(paste(c(unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Boundary"])), "comp_offset"), 
                                                 paste(c(unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Boundary"])), "comp_match_F"),
                                                 paste(c(unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Boundary"])), "comp_match_R")), 
                               reaction = rep(paste((unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Boundary"])),"composition"), times = 3), 
                               direction = c(rep("F", length(biomassS[1,])*2), rep("R", length(biomassS[1,]))))
  
  biomassConv <- list()
  for(a_rxn in biomassRxSplit$rxDesignation){
    biomassConv[[a_rxn]]$varCategory <- strsplit(a_rxn, split = ' comp')[[1]][1]
    category_species <- treatment_par[[1]][["boundaryFlux"]][[biomassConv[[a_rxn]]$varCategory]]$exchange$AltName
    category_species <- sapply(category_species, function(x){comp_met[comp_met$SpeciesName == x,]$SpeciesID[1]}) 
    met_row <- sapply(unname(category_species), function(x){c(1:length(metabolites))[metabolites == x] })
    biomassConv[[a_rxn]]$conversion <- data.frame(name = names(category_species), ID = unname(category_species), index = unname(met_row))
    
  }
  
  biomassS_split <- cbind(biomassS, biomassS, biomassS)
  colnames(biomassS_split) <- biomassRxSplit$rxDesignation
  
  S <- cbind(S_rxns_split, freeS_split, nutrientS_split, effluxS_split, biomassS_split)
  
  Sinfo <- rbind(stoiRxSplit, freeRxSplit, nutrientRxSplit, effluxRxSplit, biomassRxSplit)
  
  S <- S * t(t(rep(1, times = nrow(S)))) %*% t(ifelse(Sinfo$direction == "R", -1, 1)) #invert stoichiometry for backwards flux
  
  assign("nutrientRxSplit", nutrientRxSplit, envir=globalenv())
  assign("biomassConv", biomassConv, envir=globalenv())
  assign("Sinfo", Sinfo, envir=globalenv())
  assign("S", S, envir=globalenv())
  
}

setup_inequality_constraints <- function(){
  
  ############ Gv >= h : bounds ########
  
  ## previous splitting of reversible reactions, allows restriction of each reaction's flux to be either non-negative or non-positive (depending on constrained direction) ##
  
  ## since for the gurobi setup nutrient fluxes were split into 3 fluxes in order to accomidate the QP without an offset setup, the linear combination of uptake fluxes should instead be bounded
  ## this is difficult to do with this solver so instead for each molecule of a nutrient that is taken up an additional 'bookkeeping nutrient' will also be taken up.  The efflux of this bookkeeping
  ## nutrient through a single reaction can then be bounded by the empirical maximum nutrient influx rate, maintenance of global flux balance will enforce that the sum of actual nutrient influxes can't be higher
  ## than that of the directly bounded bookkeeping nutrient efflux
  
  trackingS <- matrix(0, ncol =  length(S[1,]), nrow = length(unique(nutrientRxSplit$reaction)))
  rownames(trackingS) <- paste(unique(nutrientRxSplit$reaction), "_bookkeeping", sep = "")
  
  for(arxn in unique(nutrientRxSplit$reaction)){
    trackingS[rownames(trackingS) == paste(arxn, "_bookkeeping", sep = ""), Sinfo$reaction == arxn] <- ifelse(Sinfo$direction[Sinfo$reaction == arxn] == "F", 1, -1)
  }
  
  ## add in a reaction which siphons off each bookkeeping nutrient
  
  S <- rbind(S, trackingS)
  
  bookkeepingRx <- data.frame(rxDesignation = paste(unique(nutrientRxSplit$reaction), "_bookkeeping", sep = ""),
                              reaction = paste(unique(nutrientRxSplit$reaction), "_bookkeeping", sep = ""), direction = "F")
  
  bookkeepingS <- matrix(0, ncol = length(unique(nutrientRxSplit$reaction)), nrow = length(S[,1]))
  colnames(bookkeepingS) <- bookkeepingRx$rxDesignation
  
  for(arxn in unique(nutrientRxSplit$reaction)){
    bookkeepingS[rownames(S) == paste(arxn, "_bookkeeping", sep = ""),colnames(bookkeepingS) == paste(arxn, "_bookkeeping", sep = "")] <- -1
  }
  
  Sinfo <- rbind(Sinfo, bookkeepingRx)
  S <- cbind(S, bookkeepingS)
  
  assign("Sinfo", Sinfo, envir=globalenv())
  assign("S", S, envir=globalenv())
  
}


#### Inferring fluxes

calculated_QP_fluxes <- function(){
  
  # Carryout quadratic programming to optimize fluxes s.t.
  # minimizing deviations between experimental boundary fluxes and flux-balanced fluxes
  
  require(Matrix) #for sparse matrix class
  require(gurobi) #interface for gurobi solver (license required: free for academics)
  #install.packages('/Library/gurobi600/mac64/R/gurobi_6.0-0.tgz', repos=NULL)
  
  qpModel <- list()
  qpparams <- list(OptimalityTol = 10^-9, FeasibilityTol = 10^-9, BarConvTol = 10^-16)
  
  flux_elevation_factor <- 10^6
  flux_penalty <- 500/(flux_elevation_factor)
  
  qpModel$A <- Matrix(S)
  qpModel$rhs <- Fzero #flux balance
  qpModel$sense <- rep("=", times = length(S[,1])) #global flux balance
  qpModel$lb <- rep(0, times = length(S[1,])) #all fluxes are greater than zero
  qpModel$ub <- rep(Inf , times = length(S[1,])) #overwrite bookkeeping fluxes with empirical maximum rates
  
  qpModel$Q <- diag(rep(0, length(S[1,]))) #min t(v)Qv
  
  qpModel$obj <- rep(flux_penalty, length(S[1,])) #min c * v where v >= 0 to 
  qpModel$obj[grep('^r_', Sinfo$rxDesignation, invert = TRUE)] <- 0
  # apply L1 penalties to corresponding reactions
  qpModel$obj[qpModel$obj != 0] <- qpModel$obj[qpModel$obj != 0] * prot_penalty[data.table::chmatch(Sinfo$reaction[qpModel$obj != 0], names(prot_penalty))]
  qpModel$obj[qpModel$obj == 0] <- 10^-8 # reduce zero L1 penalties to a minute value, to aid in finding feasible solutions
  
  ### QP-specific output_files ###
  
  flux_vectors <- list()
  
  ### iterate through conditions and optimize the fit of fluxes to boundary and biomass conditions ###
  
  for(treatment in 1:nrow(chemostatInfo)){
    
    # Define the maximal flux through reactions which are intrinsically constrained by the concentration of nutrients in the media, or are auxotrophies
    
    cond_bound <- rep(Inf, times = length(S[1,]))
    cond_bound[Sinfo$reaction %in% treatment_par[[treatment]]$auxotrophies] <- 0 #auxotrophies have a maximum flux of zero
    
    cond_nutrients <- treatment_par[[treatment]]$nutrients[treatment_par[[treatment]]$nutrients$type == "uptake",]
    cond_nutrients$index <- sapply(paste(cond_nutrients$specie, "boundary_bookkeeping"), function(x){c(1:length(Sinfo[,1]))[Sinfo$rxDesignation == x]})
    
    cond_bound[cond_nutrients$index] <- cond_nutrients$ub #maximal nutrient fluxes set as [nutrient]*DR
    
    qpModel$ub <- cond_bound #hard bound the maximal flux through each reaction - Inf except for nutrient absorption
    
    # Define lower bounds on flux through reactions - zero unless otherwise informed
    
    cond_min <- rep(0, times = length(S[1,]))
    cond_min[cond_nutrients$index] <- cond_nutrients$lb
    
    qpModel$lb <- cond_min
    
    # Setup quadratic penalty
    
    ## constrain the offset fluxes to exactly equal the expected flux ##
    # exchange rates for nutrient/excreted mets, lb = ub = Glucose uptake = fluxComp for composition
    
    cond_boundary_rxns <- data.frame(Sinfo[grep('offset', Sinfo[,1]),], index = grep('offset', Sinfo[,1]))
    cond_boundary_rxns$rate <- NA
    for(nutrient in treatment_par[[treatment]]$nutrients$specie){
      cond_boundary_rxns$rate[cond_boundary_rxns$reaction == paste(nutrient, "boundary")] <- treatment_par[[treatment]]$nutrients$change[treatment_par[[treatment]]$nutrients$specie == nutrient]
    }
    fluxComp = cond_boundary_rxns$rate[ cond_boundary_rxns$rxDesignation == 'D-glucose [extracellular] boundary_offset']
    cond_boundary_rxns$rate[cond_boundary_rxns$reaction %in% paste(unique(comp_by_cond$compositionFile$varCategory), "composition")] <- fluxComp
    
    qpModel$lb[cond_boundary_rxns$index][!is.na(cond_boundary_rxns$rate)] <- cond_boundary_rxns$rate[!is.na(cond_boundary_rxns$rate)]
    qpModel$ub[cond_boundary_rxns$index][!is.na(cond_boundary_rxns$rate)] <- cond_boundary_rxns$rate[!is.na(cond_boundary_rxns$rate)]
    
    qpModel$lb[cond_boundary_rxns$index][is.na(cond_boundary_rxns$rate)] <- 0
    qpModel$ub[cond_boundary_rxns$index][is.na(cond_boundary_rxns$rate)] <- Inf
    
    ## quadratic matching of exchange fluxes and production of biomass components ##
    
    matchedSpecies <- data.frame(Sinfo[grep('match', Sinfo[,1]),], index = grep('match', Sinfo[,1]))
    matchedSpecies$Precision <- NA
    
    ## input the precision of each media specie - uptake and excretion ##
    
    for(nutrientSpec in treatment_par[[treatment]]$nutrients$specie){
      matchedSpecies$Precision[matchedSpecies$reaction == paste(nutrientSpec, "boundary")] <-
        (1/treatment_par[[treatment]]$nutrients$sd[treatment_par[[treatment]]$nutrients$specie == nutrientSpec])^2
    }
    
    stacked_comp_offset <- rbind(matchedSpecies[,1:4], cond_boundary_rxns[,1:4])
    
    # normalize all mass reactions to fluxComp (glucose uptake)
    for(biomassSpec in names(treatment_par[[treatment]]$boundaryFlux)){
      ## overwrite diagonal elements of Q with precision (inverse variance) ##
      matchedSpecies$Precision[matchedSpecies$reaction == paste(biomassSpec, "composition")] <- (1/(treatment_par[[treatment]]$boundaryFlux[[biomassSpec]]$SD*fluxComp))^2
      biomassSpec_regex <- gsub('\\[', '\\\\[', biomassSpec)
      biomassSpec_regex <- gsub('\\]', '\\\\]', biomassSpec_regex)
      
      if(biomassSpec %in% unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Boundary"])){
        ### overwrite stoichiometry for each biomass reaction to reflect the actual portions consumed such that the expected flux is equal to fluxComp ##
        
        for(biomassSpecRx in names(biomassConv)[grep(paste0('^', biomassSpec_regex), names(biomassConv))]){
          
          comp_replacement <- data.frame(treatment_par[[treatment]]$boundaryFlux[[biomassSpec]]$exchange, biomassConv[[biomassSpecRx]])
          
          qpModel$A[comp_replacement$conversion.index,stacked_comp_offset$index[stacked_comp_offset$rxDesignation == biomassSpecRx]] <- comp_replacement$change*ifelse(sum(grep('_R$', biomassSpecRx)) != 0, -1, 1) / fluxComp
        }
      }
      
      if(biomassSpec %in% unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Internal"])){
        ### overwrite internal fluxes ##
        for(biomassSpecRx in names(biomassConv)[grep(paste('^', biomassSpec_regex, sep = ""), names(biomassConv))]){
          
          comp_replacement <- internalRxSplit[internalRxSplit$reaction %in% paste(biomassSpec, "composition"),]
          comp_replacement$change <- treatment_par[[treatment]]$boundaryFlux[[biomassSpec]]$exchange$change * ifelse(comp_replacement$direction == "F", 1, -1) / fluxComp
          
          qpModel$A[rownames(S) == paste(biomassSpec, "usage_bookkeeping"), chmatch(internalRxSplit$rxDesignation[internalRxSplit$reaction %in% paste(biomassSpec, "composition")], colnames(S))] <- comp_replacement$change
          
        }
      }
    }
    
    ## correct precisions for flux elevation factor
    matchedSpecies$Precision <- matchedSpecies$Precision / flux_elevation_factor^2
    
    ## for some species, no precision estimate is provided because they weren't directly measured.
    ## By enforcing some arbitrary penalty on these reactions the unbounded stoichiometrically equivalent 'offset fluxes' should carry flux
    matchedSpecies$Precision[is.na(matchedSpecies$Precision)] <- 1
    #matchedSpecies$Precision[is.na(matchedSpecies$Precision)] <- median(matchedSpecies$Precision,na.rm=T)
    
    diag(qpModel$Q)[matchedSpecies$index] <- matchedSpecies$Precision
    
    qpModel$lb <- qpModel$lb * flux_elevation_factor
    
    # Set the upper bounds to maximally maxMult times the glucose uptake rate
    maxMult = 100
    qpModel$ub <- sapply(qpModel$ub,function(x){min(maxMult*fluxComp, x)})
    qpModel$ub <- qpModel$ub * flux_elevation_factor
    
    
    ### Run standard quadratic programming ###
    
    solvedModel <- gurobi(qpModel, qpparams) #solve with barrier algorithm
    
    if(solvedModel$status == "NUMERIC"){
      alt_parms <- qpparams; alt_parms$method <- 1
      solvedModel <- gurobi(qpModel, alt_parms) #solve with dual simplex
    }
    
    collapsedFlux <- sapply(unique(Sinfo$reaction), function(frcombo){
      frindices <- c(1:length(Sinfo[,1]))[Sinfo$reaction == frcombo]
      sum(ifelse(Sinfo$direction[frindices] == "F", 1, -1) * solvedModel$x[frindices]/flux_elevation_factor)
    })
    
    # set fluxes below 10^-10 to zero, hist(log10(collapsedFlux)) indicates that there are many fluxes which are non-zero because of failure to round to zero
    collapsedFlux[abs(collapsedFlux) < 10^-10] <- 0
    
    # adjust composition fluxes to remove glucose uptake rate and restore units
    
    collapsedFlux[grep('composition', names(collapsedFlux))] <- sapply(grep('composition', names(collapsedFlux)), function(x){
      x <- collapsedFlux[x]
      x/fluxComp * abs(treatment_par[[names(treatment_par)[treatment]]]$"boundaryFlux"[[sub(' composition', '', names(x))]]$exchange$change[1])
    })
    
    # save summaries of flux carried and deviations between allowable and experimentally-measured fluxes
    flux_vectors[[names(treatment_par)[treatment]]]$"flux" <- collapsedFlux
  }
  
  ### Aggregate condition level flux estimates into a table
  
  rxNames <- data.frame(reactionID = unique(Sinfo$reaction), Name = unique(Sinfo$reaction))
  rxNames$Name[grep('r_[0-9]+', rxNames$Name)] = unname(rxnIDtoEnz(rxNames$Name[grep('r_[0-9]+', rxNames$Name)]))
  
  fluxMat <- matrix(NA, ncol = nrow(chemostatInfo),
                    nrow = length(flux_vectors[[1]]$flux));
  colnames(fluxMat) <- names(flux_vectors)
  
  for(i in 1:nrow(chemostatInfo)){
    fluxMat[,i] <- flux_vectors[[i]]$flux
  }
  rxNames <- rxNames[rowSums(fluxMat) != 0,]
  fluxMat <- fluxMat[rowSums(fluxMat) != 0,]
  
  fluxMat_per_cellVol <- fluxMat / t(t(rep(1, length(fluxMat[,1])))) %*% chemostatInfo$VolFrac_mean # moles per h*mL cell volume
  
  return(cbind(rxNames, data.frame(fluxMat_per_cellVol)))
  
}

read_write_python_FVA <- function(status = "write", directory = NULL){
  
  # Setup flux files so that FVA can be performed
  # Gurobi cannot currently do quadratic-constrained linear programming using R so python is used
  
  # if status == "write"
  # generate pythonDat for running qp_fba_clust.py in python (ideally on a cluster)
  # one file will be generated for each condition to allow for parallelization
  # if status == "read"
  # read pythonout files generated by qp_fba_clust.py
  
  stopifnot(status %in% c("read", "write"),
            class(directory) %in% c("NULL", "character"))
  
  ### The model for FVA is nearly identical to QP
  
  require(Matrix) #for sparse matrix class
  require(gurobi) #interface for gurobi solver (license required: free for academics)
  
  if(status == "read"){
    
    fva_summary <- NULL
    
    for(treatment in 1:nrow(chemostatInfo)){
      
      pythout = FVA_setup(qpModel, useCluster)
      
      fva_sum <- modelOut[,grep('FVA[min|max]|asID|offset|runStatus|violation', colnames(modelOut))]
      
      ## add offset back in
      fva_sum[,grep('^FVA[min|max]', colnames(fva_sum))] <- fva_sum[,grep('^FVA[min|max]', colnames(fva_sum))] + t(t(ifelse(!is.na(fva_sum$offset), fva_sum$offset, 0))) %*% rep(1, length(grep('^FVA[min|max]', colnames(fva_sum))))
      ## remove flux elevation factor
      fva_sum[,grep('^FVA[min|max]', colnames(fva_sum))] <- fva_sum[,grep('^FVA[min|max]', colnames(fva_sum))]/flux_elevation_factor
      fva_sum <- fva_sum[,colnames(fva_sum) != "offset"]
      fva_sum <- data.table(melt(fva_sum, id.vars = "asID"))
      
      fva_sum[,dataType := 
                if(length(grep('^FVA[min|max]', variable)) != 0){"optim"}else if(length(grep('^runStatus', variable)) != 0){"status"}else if(length(grep('^violation', variable)) != 0){"violation"}
              , by = variable]
      
      levels(fva_sum$variable) <- sub('[A-Za-z.]+FVA', 'FVA', levels(fva_sum$variable))
      
      fva_sum <- data.table(dcast(fva_sum, "asID + variable ~ dataType", value.var = "value"))
      fva_sum$optim <- as.numeric(fva_sum$optim)
      
      fva_sum[, dataType:= strsplit(as.character(variable), split = "_")[[1]][1], by = variable]
      fva_sum[, boundType:= strsplit(as.character(variable), split = "_")[[1]][1], by = variable]
      fva_sum[, logLikelihood:= strsplit(as.character(variable), split = "_")[[1]][2], by = variable]
      fva_sum[, relLikelihood:= round(as.numeric(logLikelihood) - max(as.numeric(logLikelihood)), 5), by = logLikelihood]
      fva_sum$treatment = treatment
      fva_sum = fva_sum[,list(treatment, asID, boundType, relLikelihood, optim, status, violation),]
      
      fva_sum$optim[grep('comp', fva_sum$asID)] <- sapply(grep('comp', fva_sum$asID), function(x){
        x = data.frame(ID = sub(' comp', '', fva_sum$asID[x]), flux = fva_sum$optim[x])
        x[1,2]/fluxComp * abs(treatment_par[[names(treatment_par)[treatment]]]$"boundaryFlux"[[x[1,1]]]$exchange$change[1])
      })
      
      fva_summary <- rbind(fva_summary, fva_sum)
    }
  }
  
  if(status == "write"){
    
    qpModel <- list()
    qpparams <- list(OptimalityTol = 10^-9, FeasibilityTol = 10^-9, BarConvTol = 10^-16)
    
    flux_elevation_factor <- 10^6
    flux_penalty <- 500/(flux_elevation_factor)
    
    qpModel$A <- Matrix(S)
    qpModel$rhs <- Fzero #flux balance
    qpModel$sense <- rep("=", times = length(S[,1])) #global flux balance
    qpModel$lb <- rep(0, times = length(S[1,])) #all fluxes are greater than zero
    qpModel$ub <- rep(Inf , times = length(S[1,])) #overwrite bookkeeping fluxes with empirical maximum rates
    
    qpModel$Q <- diag(rep(0, length(S[1,]))) #min t(v)Qv
    
    qpModel$obj <- rep(flux_penalty, length(S[1,])) #min c * v where v >= 0 to 
    qpModel$obj[grep('^r_', Sinfo$rxDesignation, invert = TRUE)] <- 0
    # apply L1 penalties to corresponding reactions
    qpModel$obj[qpModel$obj != 0] <- qpModel$obj[qpModel$obj != 0] * prot_penalty[data.table::chmatch(Sinfo$reaction[qpModel$obj != 0], names(prot_penalty))]
    qpModel$obj[qpModel$obj == 0] <- 10^-8 # reduce zero L1 penalties to a minute value, to aid in finding feasible solutions
    
    ### QP-specific output_files ###
    
    flux_vectors <- list()
    
    ### iterate through conditions and optimize the fit of fluxes to boundary and biomass conditions ###
    
    for(treatment in 1:nrow(chemostatInfo)){
      
      # Define the maximal flux through reactions which are intrinsically constrained by the concentration of nutrients in the media, or are auxotrophies
      
      cond_bound <- rep(Inf, times = length(S[1,]))
      cond_bound[Sinfo$reaction %in% treatment_par[[treatment]]$auxotrophies] <- 0 #auxotrophies have a maximum flux of zero
      
      cond_nutrients <- treatment_par[[treatment]]$nutrients[treatment_par[[treatment]]$nutrients$type == "uptake",]
      cond_nutrients$index <- sapply(paste(cond_nutrients$specie, "boundary_bookkeeping"), function(x){c(1:length(Sinfo[,1]))[Sinfo$rxDesignation == x]})
      
      cond_bound[cond_nutrients$index] <- cond_nutrients$ub #maximal nutrient fluxes set as [nutrient]*DR
      
      qpModel$ub <- cond_bound #hard bound the maximal flux through each reaction - Inf except for nutrient absorption
      
      # Define lower bounds on flux through reactions - zero unless otherwise informed
      
      cond_min <- rep(0, times = length(S[1,]))
      cond_min[cond_nutrients$index] <- cond_nutrients$lb
      
      qpModel$lb <- cond_min
      
      # Setup quadratic penalty
      
      ## constrain the offset fluxes to exactly equal the expected flux ##
      # exchange rates for nutrient/excreted mets, lb = ub = Glucose uptake = fluxComp for composition
      
      cond_boundary_rxns <- data.frame(Sinfo[grep('offset', Sinfo[,1]),], index = grep('offset', Sinfo[,1]))
      cond_boundary_rxns$rate <- NA
      for(nutrient in treatment_par[[treatment]]$nutrients$specie){
        cond_boundary_rxns$rate[cond_boundary_rxns$reaction == paste(nutrient, "boundary")] <- treatment_par[[treatment]]$nutrients$change[treatment_par[[treatment]]$nutrients$specie == nutrient]
      }
      fluxComp = cond_boundary_rxns$rate[ cond_boundary_rxns$rxDesignation == 'D-glucose [extracellular] boundary_offset']
      cond_boundary_rxns$rate[cond_boundary_rxns$reaction %in% paste(unique(comp_by_cond$compositionFile$varCategory), "composition")] <- fluxComp
      
      qpModel$lb[cond_boundary_rxns$index][!is.na(cond_boundary_rxns$rate)] <- cond_boundary_rxns$rate[!is.na(cond_boundary_rxns$rate)]
      qpModel$ub[cond_boundary_rxns$index][!is.na(cond_boundary_rxns$rate)] <- cond_boundary_rxns$rate[!is.na(cond_boundary_rxns$rate)]
      
      qpModel$lb[cond_boundary_rxns$index][is.na(cond_boundary_rxns$rate)] <- 0
      qpModel$ub[cond_boundary_rxns$index][is.na(cond_boundary_rxns$rate)] <- Inf
      
      ## quadratic matching of exchange fluxes and production of biomass components ##
      
      matchedSpecies <- data.frame(Sinfo[grep('match', Sinfo[,1]),], index = grep('match', Sinfo[,1]))
      matchedSpecies$Precision <- NA
      
      ## input the precision of each media specie - uptake and excretion ##
      
      for(nutrientSpec in treatment_par[[treatment]]$nutrients$specie){
        matchedSpecies$Precision[matchedSpecies$reaction == paste(nutrientSpec, "boundary")] <-
          (1/treatment_par[[treatment]]$nutrients$sd[treatment_par[[treatment]]$nutrients$specie == nutrientSpec])^2
      }
      
      stacked_comp_offset <- rbind(matchedSpecies[,1:4], cond_boundary_rxns[,1:4])
      
      # normalize all mass reactions to fluxComp (glucose uptake)
      for(biomassSpec in names(treatment_par[[treatment]]$boundaryFlux)){
        ## overwrite diagonal elements of Q with precision (inverse variance) ##
        matchedSpecies$Precision[matchedSpecies$reaction == paste(biomassSpec, "composition")] <- (1/(treatment_par[[treatment]]$boundaryFlux[[biomassSpec]]$SD*fluxComp))^2
        biomassSpec_regex <- gsub('\\[', '\\\\[', biomassSpec)
        biomassSpec_regex <- gsub('\\]', '\\\\]', biomassSpec_regex)
        
        if(biomassSpec %in% unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Boundary"])){
          ### overwrite stoichiometry for each biomass reaction to reflect the actual portions consumed such that the expected flux is equal to fluxComp ##
          
          for(biomassSpecRx in names(biomassConv)[grep(paste0('^', biomassSpec_regex), names(biomassConv))]){
            
            comp_replacement <- data.frame(treatment_par[[treatment]]$boundaryFlux[[biomassSpec]]$exchange, biomassConv[[biomassSpecRx]])
            
            qpModel$A[comp_replacement$conversion.index,stacked_comp_offset$index[stacked_comp_offset$rxDesignation == biomassSpecRx]] <- comp_replacement$change*ifelse(sum(grep('_R$', biomassSpecRx)) != 0, -1, 1) / fluxComp
          }
        }
        
        if(biomassSpec %in% unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Internal"])){
          ### overwrite internal fluxes ##
          for(biomassSpecRx in names(biomassConv)[grep(paste('^', biomassSpec_regex, sep = ""), names(biomassConv))]){
            
            comp_replacement <- internalRxSplit[internalRxSplit$reaction %in% paste(biomassSpec, "composition"),]
            comp_replacement$change <- treatment_par[[treatment]]$boundaryFlux[[biomassSpec]]$exchange$change * ifelse(comp_replacement$direction == "F", 1, -1) / fluxComp
            
            qpModel$A[rownames(S) == paste(biomassSpec, "usage_bookkeeping"), chmatch(internalRxSplit$rxDesignation[internalRxSplit$reaction %in% paste(biomassSpec, "composition")], colnames(S))] <- comp_replacement$change
            
          }
        }
      }
      
      ## correct precisions for flux elevation factor
      matchedSpecies$Precision <- matchedSpecies$Precision / flux_elevation_factor^2
      
      ## for some species, no precision estimate is provided because they weren't directly measured.
      ## By enforcing some arbitrary penalty on these reactions the unbounded stoichiometrically equivalent 'offset fluxes' should carry flux
      matchedSpecies$Precision[is.na(matchedSpecies$Precision)] <- 1
      #matchedSpecies$Precision[is.na(matchedSpecies$Precision)] <- median(matchedSpecies$Precision,na.rm=T)
      
      diag(qpModel$Q)[matchedSpecies$index] <- matchedSpecies$Precision
      
      qpModel$lb <- qpModel$lb * flux_elevation_factor
      
      # Set the upper bounds to maximally maxMult times the glucose uptake rate
      maxMult = 100
      qpModel$ub <- sapply(qpModel$ub,function(x){min(maxMult*fluxComp, x)})
      qpModel$ub <- qpModel$ub * flux_elevation_factor
      
      
      ### Now that the model is setup it can be saved for python
      
      FVA_setup(qpModel, directory)
      
    }
  }
}
  
FVA_setup <- function(aModel, directory){
    
    ### Setup a model that can be used to run Gurobi through python - specifically implementing flux variability analysis ###
    
    pythonMode = 'simple'
    FVA <- T
    
    # This is a character vector with length:
    # tables are sent row by row + 1 row for colnames, 1 for rownames
    # vectors as 1 row
    # each row has the variable name as first entry
    # variables: A(tab), rhs(vec), sense(vec), lb(vec), ub(vec), Q(tab), obj(vec), parameters(vec)
    # namesrxndG(vec),CI95rxndG(vec),namesthdynMet(vec),lbthdynMet(vec),ubthdynMet(vec),pythonMode(vec)  
    
    # 10 vec + 4 row/col names = 14 single rows
    lenPyt <- 15 + nrow(aModel$A) + nrow(aModel$Q)
    
    pythonDat <- vector(mode="character",length=lenPyt)
    pos =1
    
    pythonDat[pos] = paste(c('mode',pythonMode),collapse='\t')
    pos <- pos+1
    pythonDat[pos] = paste(c('FVA',FVA),collapse='\t')
    pos <- pos+1
    
    for(i in 1:nrow(aModel$A)){
      pythonDat[pos] = paste(c('A',aModel$A[i,]),collapse='\t')
      pos <- pos+1
    }
    pythonDat[pos] = paste(c('rownamesA',rownames(aModel$A)),collapse='\t')
    pos <- pos+1
    pythonDat[pos] = paste(c('colnamesA',colnames(aModel$A)),collapse='\t')
    pos <- pos+1
    
    for(i in 1:nrow(aModel$Q)){
      pythonDat[pos] = paste(c('Q',aModel$Q[i,]),collapse='\t')
      pos <- pos+1
    }
    pythonDat[pos] = paste(c('rownamesQ',rownames(aModel$Q)),collapse='\t')
    pos <- pos+1
    pythonDat[pos] = paste(c('colnamesQ',colnames(aModel$Q)),collapse='\t')
    pos <- pos+1
    
    pythonDat[pos] = paste(c('rhs',aModel$rhs),collapse='\t')
    pos <- pos+1
    pythonDat[pos] = paste(c('sense',aModel$sense),collapse='\t')
    pos <- pos+1
    pythonDat[pos] = paste(c('lb',aModel$lb),collapse='\t')
    pos <- pos+1
    pythonDat[pos] = paste(c('ub',aModel$ub),collapse='\t')
    pos <- pos+1
    pythonDat[pos] = paste(c('obj',aModel$obj),collapse='\t')
    
    write(pythonDat,file=paste0(directory, '/pythonDat_',treatment,'.txt'),sep='\t')
    
  }

FVA_read <- function(pythout, directory){
  
  pythout = (read.table(paste0(directory, '/pythonDat_',treatment,'.txt'),sep='$'))[,1]
  
  ### Parse flux variability model which was run in python ###
  
  idxS = which(pythout == 'output_start')+1
  idxE = which(pythout == 'output_end')-1
  if (length(idxE)==0){
    idxE = length(pythout)
  }
  
  solvedModel = list()
  modelOut <- data.frame(param = sapply(pythout[idxS:idxE],function(x){strsplit(x,'\t')[[1]][1]}))
  modelOut$value <- as.numeric(sapply(pythout[idxS:idxE],function(x){strsplit(x,'\t')[[1]][2]}))
  modelOut$runStatus <- sapply(pythout[idxS:idxE],function(x){sub('status_', '', strsplit(x,'\t')[[1]][3])}) # extract optimization status codes for FVA
  modelOut$violation <- as.numeric(sapply(pythout[idxS:idxE],function(x){sub('viol_', '', strsplit(x,'\t')[[1]][4])}))
  
  solvedModel$x <- modelOut$value[ modelOut$param %in% colnames(qpModel$A)]
  print(1)
  modelOut$Type <- sapply(modelOut$param,function(x){
    if (substr(x,1,5) =='RTlnc'){
      return('conc')
    } else if (substr(x,1,3) =='dGr'){
      if(substr(x,nchar(x),nchar(x))=='R'){
        return('dGrR')
      } else {
        return('dGrF')
      }
    } else if (substr(x,nchar(x)-1,nchar(x)) == 'sw'){
      return('sw')
    } else if (substr(x,nchar(x)-5,nchar(x)) == 'offset'){
      return('offset')
    } else if(substr(x,nchar(x),nchar(x))=='F'){
      return('rxnF')
    } else if(substr(x,nchar(x),nchar(x))=='R'){
      return('rxnR')
    } else if (substr(x,nchar(x)-2,nchar(x)) == 'min'){
      return(paste('FVAmin',gsub('(.*_FVA_)||(_min)','',x),sep='_'))
    } else if (substr(x,nchar(x)-2,nchar(x)) == 'max'){
      return(paste('FVAmax',gsub('(.*_FVA_)||(_max)','',x),sep='_'))
    } else return ('rxnF')
  })
  print(1.5)
  # what is the type described by the parameter (is it associated with a reaction or metabolites (eg RTln))
  modelOut$asType <- sapply(modelOut$param,function(x){
    if(substr(x,1,5) =='RTlnc'){
      return('met')
    } else {
      return('rxn')
    }
  })
  
  modelOut$asID <- sapply(1:length(modelOut$param),function(x){
    name =modelOut$param[x]
    name = gsub('_match','',name)
    if (modelOut$asType[x] == 'rxn'){
      idx = regexec('r_...._',name)[[1]][1]
      if (idx >0){
        return(substr(name,idx,idx+5))
      } else {
        if (modelOut$Type[x] == 'sw'){
          return(substr(name,1,nchar(name)-3))
        } else if (modelOut$Type[x] == 'offset'){
          return(substr(name,1,nchar(name)-7))
        } else if (modelOut$Type[x] %in% c('dGrR','dGrF')){
          return(substr(name,5,nchar(name)-2))
        } else if (grepl('FVA',modelOut$Type[x])){
          return(strsplit(name,'_FVA')[[1]][1])
        } else {
          return(substr(name,1,nchar(name)-2))
        }
      }
    }else { # metabolites
      idx = regexec('_s_',name)[[1]][1]
      return(substr(name,idx+1,nchar(name)))
    }
  })
  
  modelOut = reshape(modelOut,idvar = c('asID','asType'),timevar='Type',drop =c('param'),direction ='wide')
  colnames(modelOut) <- gsub('value\\.','',colnames(modelOut))
  modelOut <- modelOut[,colSums(!(is.na(modelOut)|modelOut=="") ) != 0] # remove optimization status codes for fluxes which were not posed as an individual optimizaiton problem
  
  modelOut$rxnNet <- (sapply(modelOut$rxnF,function(x){max(c(x,0),na.rm=T)})
                      -sapply(modelOut$rxnR,function(x){max(c(x,0),na.rm=T)}) +
                        sapply(modelOut$offset,function(x){max(c(x,0),na.rm=T)}))/flux_elevation_factor
  modelOut$rxnNet[abs(modelOut$rxnNet) < 1e-19] = 0
  modelOut$name <- sapply(modelOut$asID,function(x){
    if (substr(x,1,2) == 'r_'){
      return(rxnIDtoEnz(x))
    } else if (substr(x,1,2) == 's_'){
      return(metIDtoSpec(x))
    } else return(x)
  })
  modelOut$comparID <- sapply(modelOut$asID,function(x){
    if (substr(x,1,2) == 'r_'){
      return(rxnFile$Compartment[ rxnFile$ReactionID == x][1])
    } else (return(NA))
  })
  
  modelOut
  }


#### Functions to convert between IDs and species

metIDtoSpec <- function(meta, includeComp = F){
	if(includeComp){
    sapply(meta, function(x){
       paste(corrFile[corrFile$SpeciesID == x, colnames(corrFile) %in% c("SpeciesName", "Compartment")] , collapse = "-")
       })
    }else{
     sapply(meta, function(x){
       corrFile$SpeciesName[corrFile$SpeciesID == x]
       })
    }
}

rxnIDtoEnz <- function(rxn){
	sapply(rxn, function(x){
		rxnFile$Reaction[rxnFile$ReactionID == x][1]
		})}

rxn_search <- function(search_string, stoiMat = named_stoi, is_rxn = TRUE, index = FALSE){
	#search by metabolite or reactant and return all reactions and nonzero metabolites.
  #stoiMat rows and columns must be named with metabolite and enzyme common names: named_stoi
	if (is_rxn == TRUE){
		colz = grep(search_string, colnames(stoiMat), fixed = TRUE)
		} else {
		met = grep(search_string, rownames(stoiMat), fixed = TRUE)
		if (length(met) == 1){
			colz = c(1:length(stoiMat[1,]))[stoiMat[met,] != 0]
			} else {
			colz = c(1:length(stoiMat[1,]))[apply(stoiMat[met,], 2, is.not.zero)]
		}}
	
	if(length(colz) == 0){
		print("no hits")
		} else {
			if(index == TRUE){
				reactions[colz]
				} else {
			
			rxns = stoiMat[,colz]
			if(is.vector(rxns)){
				c(reactions[colz], rxns[rxns != 0])
				} else {
					output <- rbind(reactions[colz], rxns[apply(rxns, 1, is.not.zero),])
					colnames(output) = colnames(stoiMat)[colz]
					output
					}}
		}
	}