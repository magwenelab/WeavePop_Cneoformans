#!/usr/bin/env bash
parallel "seqkit faidx genomes-annotations/{}/predicted_proteins.fa $1 | seqkit replace -p '($)' -r ' sample={}'"  :::: results/samples.txt > proteins/$1.fa
