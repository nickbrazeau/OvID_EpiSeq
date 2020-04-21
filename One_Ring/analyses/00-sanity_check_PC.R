#..............................................................
# Sanity check pos controls
#..............................................................
library(tidyverse)
library(seqinr)
#..............................................................
# read in data
#..............................................................
fastas <- list.files("../Seq_Analyses/genomes/", pattern = ".fasta", full.names = T)
readfastas <- function(x){
  ret <- seqinr::read.fasta(x)
  names(ret) <- gsub(".fasta", "",
                     stringr::str_split_fixed(basename(x), "_", n=2)[,2])
  return(ret)
}
fastas <- lapply(fastas, readfastas)

#..............................................................
# read in primers
#..............................................................
fwd <- seqinr::read.fasta("../Seq_Analyses/targetRefSeqs_all_smpls_keep/Po/forwardPrimer.fasta")
toupper(paste(fwd[[1]], collapse = ""))
rev <- seqinr::read.fasta("../Seq_Analyses/targetRefSeqs_all_smpls_keep/Po/reversePrimer.fasta")
toupper(paste( Biostrings::reverseComplement(Biostrings::DNAString(paste(rev[[1]], collapse = "")))))


#..............................................................
# read in beds
#..............................................................
beds <- list.files("../Seq_Analyses/targetRefSeqs_all_smpls_keep/locationsByGenome/",
                   pattern = "18s.bed", full.names = T)
beds <- lapply(beds, function(x){
  ret <- readr::read_tsv(x, col_names = F)
  colnames(ret) <- c("genome", "start", "end", "det", "len", "strand")
  ret$genome <- gsub(".fasta", "",
                     stringr::str_split_fixed(basename(x), "_", n=2)[,2])
  ret <- ret[, c("genome", "start", "end", "len")]
  return(ret)


})

#..............................................................
# Subset Fastas
#..............................................................
gns <- tibble::tibble(fasta = fastas, bed = beds)
gns$bases <- purrr::pmap(gns, function(fasta, bed){
  subfasta <- fasta[[1]][(bed$start+1):bed$end] # beds are 0-based, R is 1
  subfasta <- tibble::tibble(pos = 1:length(subfasta),
                             bases = subfasta)
  return(subfasta)
})

#..............................................................
# look at output
#..............................................................
gns$name <- purrr::map(gns$fasta, names)
seqs <- gns %>%
  dplyr::select(c("name", "bases")) %>%
  tidyr::unnest(cols = c("name", "bases"))


seqs.plotObj <- seqs %>%
  dplyr::mutate(basefct = factor(bases, levels = c("a", "c", "t", "g"),
                                 labels = c("A", "C", "T", "G"))) %>%
  ggplot() +
  geom_tile(aes(x=pos, y = name, fill = basefct)) +
  scale_fill_manual("Bases", values = c("#e41a1c", "#4daf4a", "#ffff33", "#377eb8")) + # A C G T
  theme_bw() +
  ylab("") +
  xlab("")

plot(seqs.plotObj)


#..............................................................
# write out fastas
#..............................................................
fastasout <- gns$bases
fastasout <- lapply(fastasout, function(x){
  ret <- x$bases
  ret <- paste(ret, collapse = "")
})
names(fastasout) <- gns$name
dir.create("fastas_genomes")
seqinr::write.fasta(sequences = fastasout,
                    names = names(fastasout),
                    file.out = "fastas_genomes/subsetted_ref_fastas.fa")

