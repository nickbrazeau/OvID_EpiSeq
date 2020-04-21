library(tidyverse)
library(seqinr)

dir.create("results/seekdeep", recursive = T)


#..............................................................
# get haplo results
#..............................................................
ret <- readr::read_tsv("../SeekDeep_Analyses/analysis/popClustering/Po/analysis/selectedClustersInfo.tab.txt.gz")
write.csv(ret, "results/seekdeep/full_skdp_table.csv")
# relevant bits
retsimp <- ret %>%
  dplyr::select(c("s_Sample", "h_popUID", "c_AveragedFrac"))
write.csv(retsimp, "results/seekdeep/simplified_skdp_table.csv")

# write csv with seqs
retfasta <- ret %>%
  dplyr::select(c("h_popUID", "c_Consensus")) %>%
  dplyr::filter(!duplicated(.))

write_csv(retfasta, "results/seekdeep/Seq_toBlast.csv")

# write fasta
seqinr::write.fasta(sequences = retfasta$c_Consensus, names = retfasta$h_popUID,
                    file.out = "results/seekdeep/Seq_toBlast.fa")

#..............................................................
# Get extraction results
#..............................................................
extractret <- readr::read_tsv("../SeekDeep_Analyses/analysis/reports/allExtractionStats.tab.txt") %>%
  dplyr::mutate(used_reads_tot = stringr::str_split_fixed(used, "\\(", n = 2)[,1],
                used_reads_tot = as.numeric(used_reads_tot)) %>%
  dplyr::select(c("inputName", "Total", "used_reads_tot", dplyr::everything()))

processret <- readr::read_tsv("../SeekDeep_Analyses/analysis/reports/allProcessPairsCounts.tab.txt")

# seekdeep haplotypes
skdpused <- retsimp %>%
  dplyr::select(c("s_Sample")) %>%
  dplyr::mutate(skdpused = T) %>%
  dplyr::rename(inputName = s_Sample)

extractret.simp <- extractret %>%
  dplyr::left_join(., skdpused, by = "inputName") %>%
  dplyr::select(c("inputName", "Total", "used_reads_tot", "skdpused"))

write_csv(extractret.simp, "results/seekdeep/simple_extraction_results.csv")
write_csv(processret, "results/seekdeep/process_pair_extraction_results.csv")


