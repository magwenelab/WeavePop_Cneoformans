#!/usr/bin/env bash
parallel "seqkit faidx genomes-annotations/{}/predicted_cds.fa $1 | seqkit replace -p '($)' -r ' sample={}'"  :::: results/samples.txt > cds/$1.fa
