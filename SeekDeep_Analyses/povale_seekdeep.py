#!/usr/bin/python3
###############################################################################
# Purpose:  SnakeMake File to run SeekDeep for cedar ovid project
# Authors: Nick Brazeau
#Given: FASTQ
#Return: seekdeep
###############################################################################

####### Working Directory and Project Specifics ############
workdir: '/Users/nickbrazeau/Documents/GitHub/OvID_EpiSeq/SeekDeep_Analyses/'
GENOMEDIR='/Users/nickbrazeau/Documents/GitHub/OvID_EpiSeq/genomes/'

rule all:
#    input: 'gentarg_Prep_report.txt'
    input: 'seekdeep_report.txt'
#    input: 'analysis/analysis_report.txt'

rule runanalysis:
    input: 'analysis/runAnalysis.sh'
    output: 'analysis/analysis_report.txt'
    shell: 'bash {input} --numThreads 1 ; \
    echo "seekdeep analysis complete" > output'

rule rundseekdeep:
    input: symlinks='symlinks', id='idFile.tab.txt', samplenames='sampleNames.tab.txt'
    output: 'seekdeep_report.txt'
    shell: 'SeekDeep setupTarAmpAnalysis --samples {input.samplenames} \
    --outDir analysis --inputDir {input.symlinks}/ \
    --idFile {input.id} --lenCutOffs targetRefSeqs/forSeekDeep/lenCutOffs.txt \
    --overlapStatusFnp targetRefSeqs/forSeekDeep/overlapStatuses.txt \
    --refSeqsDir targetRefSeqs/forSeekDeep/refSeqs/ \
    --extraProcessClusterCmds="--fracCutOff 0.05" \
    ; echo "seekdeep run complete" > {output}'

rule extractor:
    input: id='idFile.tab.txt'
    output: 'gentarg_Prep_report.txt'
    shell: 'SeekDeep genTargetInfoFromGenomes --primers {input.id} \
     --dout targetRefSeqs \
     --genomeDir {GENOMEDIR} --selectedGenomes KF696364_PoW18s,KF696374_PoC18s,AF145334_Pf18s,AF145336_Pm18s,AF145335_Pv18s \
     --pairedEndLength 250 --overWriteDir ; echo "elucidator setup worked" > {output}'
