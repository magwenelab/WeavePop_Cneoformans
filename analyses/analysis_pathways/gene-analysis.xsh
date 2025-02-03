import click
import os
import pandas as pd
import sys
sys.path.append('/FastData/czirion/DiversityPipeline/workflow/scripts/')
import sqlfasta

# xonsh gene-analysis.xsh --gff ../results/references/VNI/VNI.gff --gene_list crypto_pathways/crypto_cAMP_genes.tsv --sv_table ../results/dataset/files/structural_variants.tsv -sn ../results/dataset/snps/VNI_variants.tsv -o VNI_cAMP -s ../workflow/scripts/sqlfasta.py -d ../results/dataset/sequences.db

@click.command()
@click.option('--gff', '-g', help='GFF file', default=None, required=True, type=click.Path(exists=True))
@click.option('--gene_list', '-l', help='Path to the file containing the list of genes to analyze', default=None, required=True, type=click.Path(exists=True))
@click.option('--sv_table', '-sv', help='Path to the TSV file with structural variants of the lineage', default=None, required=True, type=click.Path(exists=True))
@click.option('--snps_table', '-sn', help='Path to the TSV with the variants of the lineage', default=None, required=True, type=click.Path(exists=True))
@click.option('--output', '-o', help='Name for output directory.', type=click.Path())
@click.option('--script', '-s', help='Path to the script to extract sequences from the database', default=None, required=True, type=click.Path(exists=True))
@click.option('--db', '-d', help='Path to the database', default=None, required=True, type=click.Path(exists=True))

def analyze_genes(gff, gene_list, sv_table, snps_table, output, script, db):
    click.echo(f'Using GFF file: {gff}')
    click.echo(f'Using gene list: {gene_list}')
    click.echo(f'Using SV file: {sv_table}')
    click.echo(f'Using SNPs file: {snps_table}')
    os.makedirs(output, exist_ok=True)
    temp_dir = os.path.join(output, "temp") 
    os.makedirs(temp_dir, exist_ok=True)
    temp_list = os.path.join(temp_dir, "list.txt")
    temp_gff = os.path.join(temp_dir, "temp.gff")
    temp_bed = os.path.join(temp_dir, "temp.bed")
    genes_bed = os.path.join(temp_dir, "genes.bed")
    temp_sv = os.path.join(temp_dir, "temp_sv.bed")
    temp_isec_sv = os.path.join(temp_dir, "temp_isec_sv.bed")
    output_sv = os.path.join(output, "sv_genes.tsv")
    output_snps = os.path.join(output, "snps_genes.tsv")

    print("Filtering GFF file and converting to BED")
    $(cut -f1 @(gene_list) > @(temp_list))
    $(agat_sp_filter_feature_from_keep_list.pl --gff @(gff) --keep_list @(temp_list) -o @(temp_gff))
    $(agat_convert_sp_gff2bed.pl --gff @(temp_gff) -o @(temp_bed))
    $(cut -f1,2,3,4 @(temp_bed) > @(genes_bed))

    print("Intersecting BED file with SVs")
    $(tail -n +2 @(sv_table) | cut -f1,2,3,7,11,12 > @(temp_sv))
    $(bedtools intersect -a @(temp_sv) -b @(genes_bed) -wa -wb > @(temp_isec_sv))

    if os.stat(temp_isec_sv).st_size == 0:
        print("No genes from your list are in structural variant regions")
    else:
        sv_genes = pd.read_csv(temp_isec_sv, sep='\t', header=None)
        sv_genes = sv_genes[[0, 1, 2, 3, 4, 5, 9]]
        sv_genes.columns = ['Accession', 'Start', 'End', 'Structure', 'Repeat_Category', 'Sample', 'Gene']
        sv_genes.to_csv(output_sv, sep='\t', index=False)

    print("Filtering variants in genes of interest.")
    snps = pd.read_csv(snps_table, sep='\t', header=0)
    with open(temp_list, 'r') as file:
        gene_list_contents = file.read().splitlines()
    snps = snps[snps['ID'].isin(gene_list_contents)]

    if snps.empty:
        print("No variants are present in your list of genes")
    else:
        snps.to_csv(output_snps, sep='\t', index=False)

    print("Extracting sequences of genes of interest") #Currently extracts sequences from all samples, not only from the selected lineage
    output_sequences = os.path.join(output, "sequences")
    os.makedirs(output_sequences, exist_ok=True)
    for gene in gene_list_contents:
        outfna = os.path.join(output_sequences, gene + ".fna")
        outfaa = os.path.join(output_sequences, gene + ".faa")
        python @(script) lookup --seqid @(gene + "*") --seqtype DNA @(db) @(outfna)
        python @(script) lookup --seqid @(gene + "*") --seqtype PROTEIN @(db) @(outfaa)

if __name__ == '__main__':
        analyze_genes()
