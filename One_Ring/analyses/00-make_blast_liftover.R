#########################################################################
# Purpose: Trim the 18s templates for the primers
#
# Author: Nicholas F. Brazeau
#########################################################################
library(seqinr)
library(tidyverse)

#..............................................................
# read in data
#..............................................................
fa <- seqinr::read.fasta("../genomes/panplasmo_blastdb/PanPlasmo_orig.fasta")

#..................
# decoder table
#..................
liftover_names <- function(x){
  ret <- attributes(x)$Annot
  ret <- stringr::str_split_fixed(ret, " ", n = 2)
  ret[1] <- sub(">", "", ret[1])
  ret <- tibble::tibble(name = ret[1], det = ret[2])
  return(ret)
}

fa.liftover <- lapply(fa, liftover_names) %>%
  dplyr::bind_rows()

#..................
# out
#..................
seqinr::write.fasta(fa, names = names(fa), file.out = "../genomes/panplasmo_blastdb/PanPlasmo.fasta")
readr::write_csv(fa.liftover, path = "../genomes/panplasmo_blastdb/PanPlasmo_decoder.csv")

