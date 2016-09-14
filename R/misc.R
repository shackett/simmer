read_simmer_tsv <- function(folder, file){
    read.delim(system.file("extdata", folder, file, package = "simmer"))
  }
