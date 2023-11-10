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
    #sample = "SRS404449"
    reference = Path("genomes-annotations/" + sample + "/ref.fa")
    #bamfile = snakemake.input[0]
    bamfile = Path("genomes-annotations/" + sample + "/snps.bam")
    chroms = $(grep chromosome @(reference))
    chroms_list = chroms.split("\n")
    chroms_list = list(filter(None, chroms_list))
    chroms_list = [i.split(' ',1)[0] for i in chroms_list]
    chroms_list = [i.replace('>', '') for i in chroms_list]
    out_mapq = []
    out_cov = []
    for chromosome in chroms_list:
        mapq = $(samtools stats @(bamfile) @(chromosome) | grep ^MAPQ | cut -f 2-)
        mapq = pd.Series(list(mapq.split("\n")))
        mapq = chromosome + "\t" + mapq
        mapq = mapq.str.split("\t", expand = True)
        out_mapq.append(mapq)
        cov = $(samtools stats @(bamfile) @(chromosome) | grep ^COV | cut -f 2-)
        cov = pd.Series(list(cov.split("\n")))
        cov = chromosome + "\t" + cov
        cov= cov.str.split("\t", expand = True)
        out_cov.append(cov)
    quality = pd.concat(out_mapq)
    quality = quality.dropna()
    quality.columns = ["Chromosome", "MAPQ", "Count"]
    coverage = pd.concat(out_cov)
    coverage = coverage.dropna()
    coverage.columns = ["Chromosome", "Range", "Coverage", "Count"]
    mapqfile = Path("genomes-annotations/" + sample + "/mapq.tsv")
    quality.to_csv(mapqfile, index=False)
    #quality.to_csv(snakemake.output[0], index=False)
    covfile = Path("genomes-annotations/" + sample + "/cov.tsv")
    coverage.to_csv(covfile, index=False)
    #coverage.to_csv(snakemake.output[1], index=False)

if __name__ == "__main__":
    stats() # Runs the function stats