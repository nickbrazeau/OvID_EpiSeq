#########################################################################
# Purpose: Grab blast results
#########################################################################
library(tidyverse)
library(stringr)

#..............................................................
# get kmer results
#..............................................................

ovidblastresults <- list.files(path = "~/Documents/GitHub/OvID_EpiSeq/Blast_Confirm/trimmed_results/",
                               pattern = "_magblastret.gz",
                               full.names = T)

readblast <- function(path){
  # https://ncbi.github.io/magicblast/doc/output.html
  ret <- readr::read_tsv(path, skip = 3, col_names = F)
  ret <- ret[, c(1,2,3,7,8,9,10,13,16)]
  colnames(ret) <- c("read", "RefSeqMatch", "perc_identity", "query_start", "query_end", "ref_start", "ref_end", "score", "read_length")
  ret %>%
    dplyr::mutate(
      smpl = sub("_magblastret.gz", "", basename(path))
    ) %>%
    dplyr::select(c("smpl", dplyr::everything()))
}

blastout <- lapply(ovidblastresults, readblast) %>%
  dplyr::bind_rows()

# out
dir.create("results/Blast_Confirm", recursive = T)
write_csv(blastout, "results/Blast_Confirm/blast_trimmed.csv")



#..............................................................
# Summary
#..............................................................
blastout_ovale_summary <- blastout %>%
  dplyr::filter(score >= 200) %>%
  dplyr::group_by(smpl, RefSeqMatch) %>%
  dplyr::summarise(
    spec_count = dplyr::n(),
    meanpercid = mean(perc_identity),
    sdpercid = sd(perc_identity)
  ) %>%
  dplyr::ungroup()

#..................
# liftover join
#..................
decoder <- readr::read_csv("../genomes/panplasmo_blastdb/PanPlasmo_decoder.csv") %>%
  dplyr::rename(RefSeqMatch = name)

blastout_ovale_summary <- blastout_ovale_summary %>%
  dplyr::left_join(., y = decoder, by = "RefSeqMatch")

#..............................................................
# out
#..............................................................
write_csv(blastout_ovale_summary, "results/Blast_Confirm/blast_trimmed_summary.csv")


