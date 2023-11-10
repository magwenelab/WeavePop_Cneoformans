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
reference = Path(sample + "/ref.fa")
grep chromosome @(reference)
# Trying to translate the code bellow into Xonsh

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

