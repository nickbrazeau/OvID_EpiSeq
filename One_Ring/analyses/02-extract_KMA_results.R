#########################################################################
# Purpose: Grab Kmer results
#########################################################################
library(tidyverse)
library(stringr)
library(DT)

dir.create("results/kmer/", recursive = T)

#..............................................................
# get kmer results
#..............................................................

ovidkmaresults <- list.dirs(path = "~/Documents/GitHub/OvID_EpiSeq/Kmer_Confirm/")
ovidkmaresults <- ovidkmaresults[grepl("_kma_results", ovidkmaresults)]
ovidkmaresults.paths <- paste0(ovidkmaresults, "/",
                               stringr::str_split_fixed(basename(ovidkmaresults), "_", n=2)[,1],
                               ".res")

readkma <- function(path){
  smpl <- gsub(".res", "", basename(path))
  dat <- readr::read_tsv(path)

  ret <- cbind.data.frame(smpl, dat)
  colnames(ret) <- gsub("#", "", colnames(ret))

  return(ret)


}

kmaout <- lapply(ovidkmaresults.paths, readkma) %>%
  dplyr::bind_rows()

#..............................................................
# out
#..............................................................
write_csv(kmaout, "results/kmer/kmer_kma_based_results.csv")
DT::datatable(kmaout, extensions='Buttons',
              options = list(
                searching = T,
                pageLength = 20,
                dom = 'Bfrtip',
                buttons = c('csv')))


