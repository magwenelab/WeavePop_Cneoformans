#/usr/bin/env xonsh -c

from pathlib import Path
import pandas as pd
import click

@click.command()
@click.argument("sample", type=str)

def stats(sample): # Start definition of function with the sample as argument
    """This script runs 'samtools stats' on a sample for each chromosome, 
    extracts the MAPQ and COV sections for each chromosome and 
    combines the results of all chromosomes in 
    mapq.tsv and cov.tsv, respectively.
    
    SAMPLE is the sample name.
    """
    sample = "SRS404449"
    reference = Path("genomes-annotations/" + sample + "/ref.fa")
    bamfile = snakemake.input[0]
    bamfile = Path("genomes-annotations/" + sample + "/snps.bam")

    chroms = $(grep chromosome @(reference))
    chroms_list = chroms.split("\n")
    chroms_list = list(filter(None, chroms_list))
    chroms_list = [i.split(' ',1)[0] for i in chroms_list]
    chroms_list = [i.replace('>', '') for i in chroms_list]

    for chromosome in chroms_list:
        mapq = $(samtools stats @(bamfile) @(chromosome) | grep ^MAPQ | cut -f 2-)
        mapq = pd.Series(list(mapq.split("\n")))
        mapq = chromosome + "\t" + mapq
        mapq = mapq.str.split("\t", expand = True)
        cov = $(samtools stats @(bamfile) @(chromosome) | grep ^COV | cut -f 2-)
        cov = pd.Series(list(cov.split("\n")))
        cov = chromosome + "\t" + cov
        cov= cov.str.split("\t", expand = True)

# Trying to translate the code bellow into Xonsh

#!/usr/bin/env bash



grep ">" ref.fa | cut -d" " -f1 | cut -d">" -f2 | while read line
    do 
        samtools stats snps.bam $line > $line.stats
        grep ^MAPQ $line.stats | cut -f 2- > $line.mapq
        sed -e "s/^/$line\t/" -i $line.mapq
        grep ^COV $line.stats | cut -f 2- > $line.cov
        sed -e "s/^/$line\t/" -i $line.cov
    done

cat *.cov > cov.tsv
cat *.mapq > mapq.tsv
grep ">" ref.fa | cut -d" " -f1 | cut -d">" -f2 | while read line
    do     
        rm $line.mapq $line.cov
    done