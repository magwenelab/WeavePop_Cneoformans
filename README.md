# WeavePop_Cneoformans

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18600877.svg)](https://doi.org/10.5281/zenodo.18600877)

Description of running the [WeavePop](https://github.com/magwenelab/WeavePop) in the Desjardins et al. 2017 
and Ashton et al. 2019 datasets and the analyses of the results.

## Download and processing of data and config files

The starting files were obtained as described in:
* `Crypto_Desjardins/config/input_prep.ipynb`
* `Crypto_Ashton/config/input_prep.ipynb`.

### Reference genomes

The reference genome assemblies for each lineage are:
 * VNI: Strain H99 from [FungiDB release 65](https://fungidb.org/common/downloads/release-65/CneoformansH99/).  
 * VNII: Strain VNII from [GCA_022832995.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_022832995.1/).  
 * VNBI: Strain Bt22 from [PRJNA1081442](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1081442).  
 * VNBII: Strain Bt89 from [GCA_023650575.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_023650575.1/).

Located in [here](https://github.com/magwenelab/WeavePop_Cneoformans/tree/main/Crypto_Desjardins/data/references).  
The reference genomes were aligned to the VNI reference genome to obtain assemblies that have the same orientation
and to know the homology between chromosomes as described in 
[cryptococcus_reference_genomes](https://github.com/magwenelab/cryptococcus_reference_genomes/tree/gff). 

### FASTQs

The raw Illumina sequencing reads were downloaded from the NIH Sequence Read Archive (BioProject IDs Desjardins: PRJNA382844; 
Ashton: PRJEB27222 and PRJEB5282) using the workflow [download-tools](https://github.com/magwenelab/download-tools). 
In the cases where there were multiple sequencing runs for a given BioSample, we concatenated the FASTQ files of all the paired runs into one. 
The unpaired reads were ignored. The reads were cleaned with FastP.


## Steps

1) Run the analysis workflow for the Desjardins dataset:
```
conda activate snakemake
snakemake --profile Crypto_Desjardins/config/default
```
2) Run the analysis workflow for the Ashton dataset:
```
snakemake --profile Crypto_Ashton/config/default
``` 
2) Run the join_datasets workflow:

```
snakemake --profile Crypto_Desjardins_Ashton/config/default
``` 
3) Install the Conda environment for the analyses below with `mamba env create -f analyses/envs/r_env.yaml`, then open R and run `install.packages("RRphylo")` (version tested is 3.0.2).  

4) Run the analysis in  `explore_depth.qmd` (described below).

5) Fore each dataset rename the directories `02.Dataset` to `02.Dataset_all` and rerun the workflow excluding non-haploid (`samples_to_exclude: "config/non_haploids.txt"` in the `config.yaml` files) and low quality samples (`depth_quality: filter:True`). 

6) For the joined dataset rename the directory `02.Dataset` to `02.Dataset_all` and rerun the join datasets workflow.  

7) Run the analysis (described below)
* `aneuploides.qmd` 
* `tree_merge.qmd`
* `tree_plot_cnvs.qmd`
* `summary_per_lineage.qmd`

## Analyses

The working directory of all the Quarto documents is `analyses/`.  
They are rendered with the command: `quarto render analyses/scripts/<name>.qmd`.  

All the input files used in the analyses come from the input or results of WeavePop in `Crypto_Desjardins/`, 
`Crypto_Ashton/` or `Crypto_Desjardins_Ashton/` or from the following external data.  

External data:  
* `data/raw/CryptoDiversity_Desjardins_Tree.tre`: From [CryptoDiversity_Tree_Info](https://github.com/magwenelab/CryptoDiversity_Tree_Info/blob/main/CryptoDiversity_Desjardins_Tree.tre)
* `data/raw/2017.06.09.all_ours_and_desj.snp_sites.mod.fa.cln.tree`: Received from Philip Ashton (Dec 5 2024).
<details>
<summary> Analyses files and description </summary> 

```
analyses/
├── data # Raw data and data processed by the scripts here
├── envs # `r_env.yaml` environment definition
├── notebooks # Rendered versions of Quarto documents
├── results # Results of the analyses
├── scripts # Scripts, Jupyter notebooks, and Quarto documents
```

| Analysis | Script <br /> `scripts/` | Output | Description |
|-----------------|-----------------|-----------------| -----------------|
| Explore Depth Profile of all Samples | `explore_depth.qmd` | `results/tables/ploidy.tsv`<br />`Crypto_Ashton/config/non_haploids.txt`<br />`Crypto_Desjardins/config/non_haploids.txt`  | Explore the depth plots to identify putative non-haploid samples to exclude from the analyses. |
| Tree building | `merge_trees.qmd` | `data/processed/tree_ashton.newick`<br />  `data/processed/tree_desjardins.newick`<br />  `data/processed/tree_merged.newick`<br /> Plots in `results/trees/`<br /> `data/processed/metadata_all.csv` | Merge the trees of the Ashton and Desjardins datasets. |
| Discover aneuploidies | `aneuploidies.qmd` | `results/tables/chromosome_cnv_categories.tsv`<br /> Plots in `results/figs/`| Categorize chromosomes by coverage of CNVs |
| Plot duplications in tree |`tree_plot_cnvs.qmd`| Plots in `results/trees_dups/`| Plot the merged tree with a heatmap of duplicated chromosomes.|
|Metrics of quality and variants of final dataset | `summary_per_lineage.qmd` | `data/processed/snp_counts_desjardins.csv`<br /> `data/processed/snp_counts_ashton.csv` <br /> `results/tables/per_lineage_summary_stats.tsv`| Create summary table of  mapping stats, number of variants, and CNVs. |

</details>
