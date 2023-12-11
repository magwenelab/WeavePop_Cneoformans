#/usr/bin/env xonsh

# Get MAT loci protein_coding_gene IDs: Everything between SXI1 and STE12
# awk '/SXI1/,/STE12/' references/FungiDB-65_CneoformansH99.gff.tsv | grep protein_coding_gene | cut -f13 > results/MAT.txt
# Get rRNA IDs from the level2 primary_tag = rRNA and converting the ID into gene_ID (because the level1 primary tag is ncRNA_gene and the pattern rRNA is also in descriptions that are nor rRNA genes)
# awk '$3 == "rRNA" {print $0}' references/FungiDB-65_CneoformansH99.gff.tsv  | cut -f13 | cut -d'-' -f1 > results/rRNA.txt
# Get centromere delimiting-gene IDs from Janbon 2014 

import pandas as pd
from pathlib import Path
import click

@click.command()
@click.option('--genefile', '-g', multiple=True, help = 'Path to text file with list of gene IDs (each in new line) of the locus of interest, file name will be used as locus name.')
@click.option('--referencetsv', '-r', multiple=True, help = 'Path to TSV annotation file of reference genome.')

def getloci(genefile, referencetsv):
    """This script creates an annotation table "loci_interest.tsv" with the columns: 
    "seq_id", "primary_tag", "start", "description", "gene_id", "ID", "Name", "Loci" for each gene in <genefile>. 
    It takes as many gene files as desired, each with the IDs of the genes in a locus of interest.
    The output table has the annotation of the genes in all the reference genomes given as input in <referencetsv>."""
    mydfs = []
    for lin in referencetsv:
        mydfs.append(pd.read_csv(Path(lin), sep='\t', header=0))
    annotations = pd.concat(mydfs)
    genes = annotations[["seq_id", "primary_tag", "start", "description", "gene_id", "ID", "Name"]]
    level1_bool = genes.primary_tag.str.contains("protein_coding_gene") | genes.primary_tag.str.contains("ncRNA_gene")
    level1 = genes[level1_bool]

    mylocus = []
    for locus in genefile:
        name = locus.split('/')[-1].split('.')[-2]
        gene_list = list(pd.read_csv(Path(locus), header = None)[0])
        gene = level1[level1['ID'].isin(gene_list)]
        gene.loc[:,'Loci'] = name
        mylocus.append(gene)
    all_loci = pd.concat(mylocus)

    all_loci.to_csv('results/loci_interest.tsv',index= False, sep ='\t')

if __name__ == "__main__":
    getloci()