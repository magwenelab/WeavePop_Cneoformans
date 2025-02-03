import pandas as pd
import glob

duplications = pd.read_csv('/FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/results/duplications.tsv', sep='\t', header =0)
duplications['sample_chromosome'] = duplications['sample'] + '_' + duplications['chromosome']

in_cnvs = duplications[duplications['percent_cnv_size'] >= 80]
not_in_cnvs = duplications[duplications['percent_cnv_size'] < 80]

$(mkdir -p /FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/results/duplicated_plots)
for sample_chromosome in in_cnvs['sample_chromosome']:
    samp_name = in_cnvs[in_cnvs['sample_chromosome'] == sample_chromosome]['sample'].values[0]
    dataset_name = in_cnvs[in_cnvs['sample_chromosome'] == sample_chromosome]['dataset'].values[0]
    chrom_name = in_cnvs[in_cnvs['sample_chromosome'] == sample_chromosome]['chromosome'].values[0]
    print(f"Sample: {samp_name}, Dataset: {dataset_name}, Chromosome: {chrom_name}")
    src_pattern = f'/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_{dataset_name}/results_*/01.Samples/plots/{samp_name}/depth_by_windows.png'
    src_files = glob.glob(src_pattern)
    if src_files:
        src = src_files[0]
        dst = f'/FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/results/duplicated_plots/{samp_name}_{chrom_name}.png'
        print(f"Copying {src} to {dst}")
        $(scp @(src) @(dst))
    else:
        print(f"No files found for pattern: {src_pattern}")

$(mkdir -p /FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/results/duplicated_plots_not_in_cnvs)
for sample_chromosome in not_in_cnvs['sample_chromosome']:
    samp_name = not_in_cnvs[not_in_cnvs['sample_chromosome'] == sample_chromosome]['sample'].values[0]
    dataset_name = not_in_cnvs[not_in_cnvs['sample_chromosome'] == sample_chromosome]['dataset'].values[0]
    chrom_name = not_in_cnvs[not_in_cnvs['sample_chromosome'] == sample_chromosome]['chromosome'].values[0]
    print(f"Sample: {samp_name}, Dataset: {dataset_name}, Chromosome: {chrom_name}")
    src_pattern = f'/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_{dataset_name}/results_*/01.Samples/plots/{samp_name}/depth_by_windows.png'
    src_files = glob.glob(src_pattern)
    if src_files:
        src = src_files[0]
        dst = f'/FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/results/duplicated_plots_not_in_cnvs/{samp_name}_{chrom_name}.png'
        print(f"Copying {src} to {dst}")
        $(scp @(src) @(dst))
    else:
        print(f"No files found for pattern: {src_pattern}")
