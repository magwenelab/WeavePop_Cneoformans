#/usr/bin/env bash
  
cut -d',' -f4 lineage_references.csv | tail -n4 | while read file
    do 
        lineage=$(echo $file | cut -d'.' -f1)
        grep ">" references/$file | cut -d',' -f1 | cut -d'>' -f2 |  while read line
        do 
            accession=$(echo $line | cut -d' ' -f1)
            name=$(echo $line | rev | cut -d' ' -f1 | rev)
            echo $lineage,$accession,$name
        done
    done > chromosome_names.csv

sed -i 's/contig/2\./g' chromosome_names.csv # To edit the names of the chromosome 2 contigs of Crypto VNBI
sed -i 's/(3;11)/3/g' chromosome_names.csv 
sed -i 's/(11;3)/11/g' chromosome_names.csv 
sed -i 's/,0/,/g' chromosome_names.csv 