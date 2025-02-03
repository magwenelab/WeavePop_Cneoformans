import pandas as pd
ashton = pd.read_csv('/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Ashton/results_241203/02.Dataset/depth_quality/mapping_stats.tsv', sep='\t', header =0)
bad_high_mapq = ashton[ashton['percent_high_mapq'] < 95]
samples_bad_high_mapq = bad_high_mapq['sample'].tolist()

$(mkdir -p /FastData/czirion/Crypto_Diversity_Pipeline/analyses/quality/bad_high_mapq)
for sample in samples_bad_high_mapq:
    src = f'/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Ashton/results_241203/01.Samples/plots/{sample}/depth_by_windows.png'
    dst = f'/FastData/czirion/Crypto_Diversity_Pipeline/analyses/quality/bad_high_mapq/{sample}.png'
    $(scp @(src) @(dst))

bad_properly_paired = ashton[ashton['percent_properly_paired'] < 50]
samples_bad_properly_paired = bad_properly_paired['sample'].tolist()

$(mkdir -p /FastData/czirion/Crypto_Diversity_Pipeline/analyses/quality/bad_properly_paired)
for sample in samples_bad_properly_paired:
    src = f'/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Ashton/results_241203/01.Samples/plots/{sample}/depth_by_windows.png'
    dst = f'/FastData/czirion/Crypto_Diversity_Pipeline/analyses/quality/bad_properly_paired/{sample}.png'
    $(scp @(src) @(dst))


chrom_depth_complete = pd.read_csv('/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Ashton/results_241203/04.Intermediate_files/02.Dataset/depth_quality/depth_by_chrom_good.tsv', sep='\t', header =0)
chrom_depth = chrom_depth_complete[['sample', 'global_median', 'global_mean', 'global_mode']]
chrom_depth.drop_duplicates(inplace=True)
chrom_depth['diff_median'] = (chrom_depth['global_median'] - chrom_depth['global_mode']).abs()
chrom_depth.sort_values(by='diff_median', ascending=False, inplace=True)

big_diff_median= chrom_depth[chrom_depth['diff_median'] > 2]

samples_big_diff_median = big_diff_median['sample'].tolist()

$(mkdir -p /FastData/czirion/Crypto_Diversity_Pipeline/analyses/quality/big_diff_median)
for sample in samples_big_diff_median:
    src = f'/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Ashton/results_241203/01.Samples/plots/{sample}/depth_by_windows.png'
    dst = f'/FastData/czirion/Crypto_Diversity_Pipeline/analyses/quality/big_diff_median/{sample}.png'
    $(scp @(src) @(dst))

mapping_stats_big_diff_median = ashton[ashton['sample'].isin(samples_big_diff_median)]

chrom_depth['diff_mean'] = (chrom_depth['global_mean'] - chrom_depth['global_mode']).abs()

big_diff_mean= chrom_depth[chrom_depth['diff_mean'] > 2]

samples_big_diff_mean = big_diff_mean['sample'].tolist()

$(mkdir -p /FastData/czirion/Crypto_Diversity_Pipeline/analyses/quality/big_diff_mean)
for sample in samples_big_diff_mean:
    src = f'/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Ashton/results_241203/01.Samples/plots/{sample}/depth_by_windows.png'
    dst = f'/FastData/czirion/Crypto_Diversity_Pipeline/analyses/quality/big_diff_mean/{sample}.png'
    $(scp @(src) @(dst))

big_both = chrom_depth[(chrom_depth['diff_mean'] > 2) & (chrom_depth['diff_median'] > 2)]
len(big_both)

#
 
duplicated = chrom_depth_complete[chrom_depth_complete['norm_chrom_median'] > 1.55]
samples_duplicated = duplicated['sample'].tolist()

$(mkdir -p /FastData/czirion/Crypto_Diversity_Pipeline/analyses/quality/duplicated)
for sample in samples_duplicated:
    src = f'/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Ashton/results_241203/01.Samples/plots/{sample}/depth_by_windows.png'
    dst = f'/FastData/czirion/Crypto_Diversity_Pipeline/analyses/quality/duplicated/{sample}.png'
    $(scp @(src) @(dst))
