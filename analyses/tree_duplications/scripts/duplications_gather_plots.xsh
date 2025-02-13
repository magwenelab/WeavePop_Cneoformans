import pandas as pd
import glob

## Get the normalized read depth by window plots for the putative duplications to analyze by hand
duplications = pd.read_csv('/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/duplications_putative.tsv', sep='\t', header =0)
duplications['sample_chromosome'] = duplications['sample'] + '_' + duplications['chromosome']

# Get the duplications that detected by called CNVs in at least 80% of the size of the chromosome
in_cnvs = duplications[duplications['percent_cnv_size'] >= 80]

$(mkdir -p /FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/duplicated_plots)
for sample_chromosome in in_cnvs['sample_chromosome']:
    samp_name = in_cnvs[in_cnvs['sample_chromosome'] == sample_chromosome]['sample'].values[0]
    dataset_name = in_cnvs[in_cnvs['sample_chromosome'] == sample_chromosome]['dataset'].values[0]
    chrom_name = in_cnvs[in_cnvs['sample_chromosome'] == sample_chromosome]['chromosome'].values[0]
    print(f"Sample: {samp_name}, Dataset: {dataset_name}, Chromosome: {chrom_name}")
    src_pattern = f'/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_{dataset_name}/results/01.Samples/plots/{samp_name}/depth_by_windows.png'
    src_files = glob.glob(src_pattern)
    if src_files:
        src = src_files[0]
        dst = f'/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/duplicated_plots/{samp_name}_{chrom_name}.png'
        print(f"Copying {src} to {dst}")
        $(scp @(src) @(dst))
    else:
        print(f"No files found for pattern: {src_pattern}")

# Get the duplications that were NOT detected by called CNVs in less than 80% of the size of the chromosome
# but by median chromosomal normalized depth only
not_in_cnvs = duplications[duplications['percent_cnv_size'] < 80]

$(mkdir -p /FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/duplicated_plots_not_in_cnvs)
for sample_chromosome in not_in_cnvs['sample_chromosome']:
    samp_name = not_in_cnvs[not_in_cnvs['sample_chromosome'] == sample_chromosome]['sample'].values[0]
    dataset_name = not_in_cnvs[not_in_cnvs['sample_chromosome'] == sample_chromosome]['dataset'].values[0]
    chrom_name = not_in_cnvs[not_in_cnvs['sample_chromosome'] == sample_chromosome]['chromosome'].values[0]
    print(f"Sample: {samp_name}, Dataset: {dataset_name}, Chromosome: {chrom_name}")
    src_pattern = f'/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_{dataset_name}/results/01.Samples/plots/{samp_name}/depth_by_windows.png'
    src_files = glob.glob(src_pattern)
    if src_files:
        src = src_files[0]
        dst = f'/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/duplicated_plots_not_in_cnvs/{samp_name}_{chrom_name}.png'
        print(f"Copying {src} to {dst}")
        $(scp @(src) @(dst))
    else:
        print(f"No files found for pattern: {src_pattern}")
