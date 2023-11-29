#!/usr/bin/env bash

mkdir cds/
cat ${snakemake_input[0]} | while read protein
do
   parallel -j ${snakemake_threads[0]} "seqkit faidx genomes-annotations/{}/predicted_cds.fa $protein | seqkit replace -p '($)' -r ' sample={}'"  :::: ${snakemake_input[1]} > cds/$protein.fa
done 2> ${snakemake_log[0]}
touch ${snakemake_output[0]}

# mkdir cds/
# cat results/protein_list.txt | while read protein
# do
#     parallel -j 4 "seqkit faidx genomes-annotations/{}/predicted_cds.fa $protein | seqkit replace -p '($)' -r ' sample={}'"  :::: samples.txt > cds/$protein.fa
# done 2> logs/cds/cds.log
# touch cds/done.txt