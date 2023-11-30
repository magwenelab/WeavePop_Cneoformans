#!/usr/bin/env bash

mkdir proteins
cat ${snakemake_input[0]} | while read protein
do
   parallel "seqkit faidx genomes-annotations/{}/predicted_proteins.fa $protein | seqkit replace -p '($)' -r ' sample={}'"  :::: ${snakemake_input[1]} > proteins/$protein.fa
done 2> ${snakemake_log[0]}
touch ${snakemake_output[0]}

# mkdir proteins/
# cat results/protein_list.txt | while read protein
# do
#     parallel -j 4 "seqkit faidx genomes-annotations/{}/predicted_proteins.fa $protein | seqkit replace -p '($)' -r ' sample={}'"  :::: samples.txt > proteins/$protein.fa
# done 2> logs/proteins/proteins.log
# touch proteins/done.txt