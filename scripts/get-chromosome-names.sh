#/usr/bin/env bash

echo "lineage,accession,chromosome" > chromosome-names.csv
name=0
for file in VNI.fasta VNBI.fasta
    do 
        lineage=$(echo $file | cut -d'.' -f1)
        grep ">" data/references/$file | cut -d',' -f1 | cut -d'>' -f2 |  while read line
        do 
            accession=$(echo $line | cut -d' ' -f1)
            name=$((name+1))
            echo $lineage,$accession,$name
        done
    done >> chromosome-names.csv

for file in  VNII.fasta VNBII.fasta
    do 
        lineage=$(echo $file | cut -d'.' -f1)
        grep ">" data/references/$file | cut -d',' -f1 | cut -d'>' -f2 |  while read line
        do 
            accession=$(echo $line | cut -d' ' -f1)
            name=$(echo $line | rev | cut -d' ' -f1 | rev)
            echo $lineage,$accession,$name
        done
    done >> chromosome-names.csv

