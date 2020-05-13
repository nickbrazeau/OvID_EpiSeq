#########################################################################
# Purpose: Trim the 18s templates for the primers
#
# Author: Nicholas F. Brazeau
#########################################################################
library(Biostrings)
library(tidyverse)

#..............................................................
# read in data
#..............................................................
fastas <- list.files("../genomes/", pattern = ".fasta", full.names = T)
fastas <- lapply(fastas, Biostrings::readDNAStringSet)
names <- lapply(fastas, function(x){
  nm <- names(x)
  nm <- stringr::str_split(nm, " ", simplify = T)
  if (any(grepl("ovale", nm)) ) {
    paste(nm[ which(nm == "Plasmodium") + 1 ], nm[ which(nm == "Plasmodium") + 2 ])
  } else {
    nm[ which(nm == "Plasmodium") + 1 ]
  }
})
names <- unlist(names)

#..............................................................
# read in primers
#..............................................................
fwd <- Biostrings::readDNAStringSet("../Kmer_Confirm/ovid_primers/ovid_forward.fa")
rev <- Biostrings::reverseComplement(Biostrings::readDNAStringSet("../Kmer_Confirm/ovid_primers/ovid_reverse.fa"))

#..............................................................
# subset fastas
#..............................................................
find_primer_pos <- function(seq, fwd, rev){
  seq <- seq[[1]]
  pos <- Biostrings::matchLRPatterns(
    Lpattern = fwd,
    Rpattern = rev,
    subject = seq,
    max.gaplength = 300)
  # out
  seq[pos@ranges@start:(pos@ranges@start+pos@ranges@width - 1)]
}


trim_primer <- function(seq, fwd, rev) {
  Biostrings::trimLRPatterns(Lpattern = fwd, Rpattern = rev, subject = seq)
}

fastas.sub <- lapply(fastas, find_primer_pos, fwd = fwd[[1]], rev = rev[[1]])
fastas.trim <- lapply(fastas.sub, trim_primer, fwd = fwd[[1]], rev = rev[[1]])


#..............................................................
# write out fastas for Kmer
#..............................................................
dir.create("../genomes/panplasmo_db/")
seqinr::write.fasta(sequences = fastas.trim,
                    names = names,
                    file.out = "../genomes/panplasmo_db/panplasmo.fa")

