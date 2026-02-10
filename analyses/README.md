## scripts/
Source files of analyses notebooks.

## notebooks/scripts/
Rendered analyses notebooks in HTML. 

## results/figures/

Figures used in WeavePop paper.

## results/figures_extra/ 

Extra figures generated in analyses

## results/tables/per_lineage_summary_stats.tsv

`lineage`: Lineage of the reference genome.  
`percent_mapped_mean`: Mean percentage of mapped reads of the samples in the corresponding lineage.  
`percent_high_mapq_mean`: Mean percentage of mapped reads with high mapping quality of the samples in the corresponding lineage.  
`coverage_good_mean`: Mean coverage of good quality reads along the total length of the reference genome of the samples in the corresponding lineage.  
`unique_variants`: Total number of small variants across the samples of the corresponding lineage.  
`n_variants_mean`: Mean number of small variants per sample of the samples in the corresponding lineage.  
`n_snps_ref_strain`: Number of small variants in the samples that belong to the same strain as the `reference genome`: Bt22 in VNBI and Bt89 in VNBII.  
`n_deletions_mean`: Mean number of deleted regions per sample of the samples in the corresponding lineage.  
`percent_coverage_deletions_mean`: Mean percentage of each sample covered by deletions in the corresponding lineage.  
`n_duplications_mean`: Mean number of duplicated regions per sample of the samples in the corresponding lineage.  
`percent_coverage_deletions_mean`: Mean percentage of each sample covered by duplications in the corresponding lineage.  
`length`: Total length in base pairs of the reference genome.  

## results/tables/chromosome_cnv_categories.tsv

Summary of metrics of CNV regions per type, chromosome and sample.

`sample`: Sample ID.  
`cnv`: Type of CNV.  
`accession`: Chromosome accession.  
`n_cnvs`: Number of CNV regions.  
`total_size`: Sum of size in bp of CNV regions.  
`coverage_percent`: Percent of length of chromosome covered by CNV regions.  
`span_percent`: Percent of length of chromosome in between leftmost and right most CNV regions in chromosome.  
`size_smallest`: Size of smallest CNV region.  
`size_largest`: Size of largest CNV region.  
`std_regions_size`: Standard deviation of size of CNV regions.  
`norm_depth_mean`: Mean of normalized depth of windows in CNV regions.  
`norm_depth_median`: Median of normalized depth of windows in CNV regions.  
`smooth_depth_mean`: Mean of smoothed normalized depth of windows in CNV regions.  
`smooth_depth_median`: Median of smoothed normalized depth of windows in CNV regions.  
`chrom_depth`: Median of depth of all windows in chromosome.  
`chrom_norm_depth`: Normalized median of depth of all windows in chromosome.  
`genome_depth`: Median of depth of all windows in genome.  
`chromosome`: Name of chromosome.  
`length`: Length of chromosome in bp.  
`lineage`: Lineage of the reference genome the sample was aligned to.  
`chrom_category_aneuploidy`: Category of chromosome, according to CNV coverage percent thresholds (Full: 80-100, Partial: 20-80, Euploid:0-20).  
`chrom_category_size`: Category of chromosome, according to CNV coverage percent thresholds (Large: 80-100, Medium: 20-80, Small:5-20, Absent: 0-5).  
`strain`: Strain name.  
`source`: Isolation source.  
`dataset`: Dataset of origin.  
`sublineage`: Sublineage (first subdivision of lineage VNI).  
`samples_in_dataset_lineage`: Number of samples in corresponding dataset and lineage.   
`samples_in_lineage`: Number of samples in corresponding lineage.  
`samples_in_sublineage`: Number of samples in corresponding sublineage.  
`samples_in_lineage_source`: Number of samples in corresponding lineage and isolation source.  
`total_samples`: Total number of samples.    

## results/tables/snp_counts_ashton.tsv and snp_counts_desjardins.tsv

Columns in order: sample ID, number of variants in `01.Samples/mapping_and_variants/<sample>/snps.raw.vcf`, number of variants in `01.Samples/mapping_and_variants/<sample>/snps.filt.vcf`.
