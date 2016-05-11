options(stringsAsFactors = F)

###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@
##### Load files describing valid reactions, species (their composition) both from the core SBML model and supplemented manual annotations #####
###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@

load_metabolic_model()

###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@
##### Load files describing boundary conditions and reaction reversibility #####
###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@

format_boundary_conditions()


format_boundary_conditions <- function(){
  
  #compositionFile <- read.delim("../Yeast_comp_energy.txt") #energy required to assimilate biomass components
  nutrientFile <- read.delim("companionFiles/flux_input_data/Boer_nutrients.txt")[1:6,1:6]; nutrientCode <- data.frame(nutrient = colnames(nutrientFile)[-1], shorthand = c("N", "P", "C", "L", "U"))
  rownames(nutrientFile) <- nutrientFile[,1]; nutrientFile <- nutrientFile[,-1]
  nutModelNames <- data.frame(commonName = rownames(nutrientFile), modelName = sapply(rownames(nutrientFile), function(x){paste(x, '[extracellular]')}))
  rownames(nutrientFile) <- nutModelNames$modelName

  load("companionFiles/flux_input_data/boundaryFluxes.Rdata") #load condition specific boundary fluxes and chemostat info (actual culture DR)

  
  
  
  
  
  }



######### Setting up stoichiometric inputs

load_metabolic_model <- function(){
  
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
  
  rxnFile <<- rbind(rxnFile, customList$rxnFile)
  rxnparFile <<- rbind(rxnparFile, customList$rxnparFile)
  corrFile <<- rbind(corrFile, customList$corrFile)
  specparFile <<- rbind(specparFile, customList$specparFile)
  fluxDirFile <<- rbind(fluxDirFile, customList$fluxDirFile)
  
  ### Determine unique metabolites and reactions ###
  
  reactions <<- sort(unique(rxnFile$ReactionID))
  rxnStoi <- rxnFile[is.na(rxnFile$StoiCoef) == FALSE,]
  metabolites <<- sort(unique(rxnStoi$Metabolite))
  
  # generate the stoichiometric matrix from 
  stoiMat <<- build_stoiMat(metabolites, reactions, corrFile, rxnFile, internal_names = TRUE)
  
  # generated a named version of the stoichiometric matrix (used?)
  #named_stoi <- stoiMat
  #met_dict <- metIDtoSpec(rownames(named_stoi)); met_dict <- sapply(c(1:length(named_stoi[,1])), function(x){met_dict[x][[1]]})
  #rxn_dict <- rxnIDtoEnz(colnames(named_stoi)); rxn_dict <- sapply(c(1:length(named_stoi[1,])), function(x){rxn_dict[x][[1]]})
  #rownames(named_stoi) <- met_dict; colnames(named_stoi) <- rxn_dict

}

parse_custom = function(customRx){
  
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

build_stoiMat = function(metabolites, reactions, corrFile, rxnFile, internal_names = FALSE){
  
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

######### Functions to convert between IDs and species ######

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
