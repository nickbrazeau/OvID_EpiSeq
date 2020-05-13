#! /bin/bash
WD=/Users/nickbrazeau/Documents/GitHub/OvID_EpiSeq/Blast_Confirm/
WAIT=30 # number of seconds to wait for files to appear, absorbing some file system latency

snakemake \
	--snakefile /Users/nickbrazeau/Documents/GitHub/OvID_EpiSeq/Blast_Confirm/magicblast_wrapper.snake \
	--configfile config_blast.yaml \
	--printshellcmds \
	--directory $WD \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT \
#	--dryrun -p \
