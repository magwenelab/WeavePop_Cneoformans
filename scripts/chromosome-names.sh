#/usr/bin/env bash
  
cut -d',' -f4 ${snakemake_input[0]} | tail -n4 | while read file
    do 
        lineage=$(echo $file | cut -d'.' -f1)
        grep ">" ${snakemake_params[refdir]}$file | cut -d',' -f1 | cut -d'>' -f2 |  while read line
        do 
            accession=$(echo $line | cut -d' ' -f1)
            name=$(echo $line | rev | cut -d' ' -f1 | rev)
            echo $lineage,$accession,$name
        done
    done > ${snakemake_output[0]} 2> ${snakemake_log[0]}

#sed -i 's/contig/2 contig/g' ${snakemake_output[0]} # To edit the names of the chromosome 2 contigs of Crypto VNBI