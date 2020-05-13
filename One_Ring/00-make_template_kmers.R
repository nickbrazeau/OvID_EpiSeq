#########################################################################
# Purpose: Trim the 18s templates for the primers
#
# Author: Nicholas F. Brazeau
#########################################################################
library(Biostrings)
library(tidyverse)

#..............................................................
# read in
#..............................................................
ref18Sseqs <- list.files("../genomes/", pattern = ".fasta", full.names = T)
ref18Sseqs <- lapply(ref18Sseqs, function(x) {Biostrings::readDNAStringSet(filepath = x, format="fasta")})

#..............................................................
# trim L and R primer
#..............................................................
Lprimer <- Biostrings::readDNAStringSet(filepath=forwardprimerpath, format = "fasta")
Rprimer <- Biostrings::reverseComplement(Biostrings::readDNAStringSet(filepath=reverseprimerpath, format = "fasta"))

Pf3D7haplotypeRef <- Biostrings::trimLRPatterns(Lpattern = Lprimer[[1]],
                                                Rpattern = Rprimer[[1]],
                                                subject = ampliconrefseq)  # trim off primers
