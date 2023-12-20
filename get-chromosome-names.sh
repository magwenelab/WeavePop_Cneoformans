#/usr/bin/env bash
  
cut -d',' -f4 files/lineage_references.csv | tail -n4 | while read file
    do 
        lineage=$(echo $file | cut -d'.' -f1)
        grep ">" references/$file | cut -d',' -f1 | cut -d'>' -f2 |  while read line
        do 
            accession=$(echo $line | cut -d' ' -f1)
            name=$(echo $line | rev | cut -d' ' -f1 | rev)
            echo $lineage,$accession,$name
        done
    done > tmp

sed -i 's/contig/2\./g' tmp # To edit the names of the chromosome 2 contigs of Crypto VNBI
sed -i 's/(3;11)/3/g' tmp
sed -i 's/(11;3)/11/g' tmp 
sed -i 's/,0/,/g' tmp 
awk -F ',' '{ gsub("0000", "0", $3) ; print }' tmp > files/chromosome_names.csv
rm tmp