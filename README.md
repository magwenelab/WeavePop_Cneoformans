# WeavePop_Cneoformans

Description of running the [WeavePop](https://github.com/magwenelab/WeavePop) in the Desjardins et al. 2017 
and Ashton et al. 2019 datasets and the analyses of the results.

## Download and processing of data and config files

### References

The reference genome assemblies for each lineage are:
 * VNI: Strain H99 from [FungiDB release 65](https://fungidb.org/common/downloads/release-65/CneoformansH99/).  
 * VNII: Strain VNII from [GCA_022832995.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_022832995.1/).  
 * VNBI: Strain Bt22 from [PRJNA1081442](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1081442).  
 * VNBII: Strain Bt89 from [GCA_023650575.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_023650575.1/).

The reference genomes were aligned to the VNI reference genome to obtain assemblies that have the same orientation
and to know the homology between chromosomes as described in 
[cryptococcus_reference_genomes](https://github.com/magwenelab/cryptococcus_reference_genomes). 
The `config/chromosomes.csv` files in both `Crypto_Desjardins/` and `Crypto_Ashton/` 
were created from the results of that pipeline with `scripts/prepare_ref_genomes.ipynb`.

The script `scripts/remove_chromosome.sh` was used to remove the mitochondrial chromosome of the `VNI.fasta` genome. 

### FASTQs

The raw Illumina sequencing reads were downloaded from the NIH Sequence Read Archive (BioProject IDs Desjardins: PRJNA382844; 
Ashton: PRJEB27222 and PRJEB5282) using the workflow [download-tools](https://github.com/magwenelab/download-tools). 
In the cases where there were multiple sequencing runs for a given BioSample, we concatenated the FASTQ files of all the paired runs into one. 
This creates the FASTQ files and the tables `config/reads_table.csv` -- Columns are `sample` (SRS ID), `run` (SRR ID), 
`file1` (read pair 1), `file2` (read pair 2), `file_unpaired` (unpaired reads). The unpaired reads were ignored.

The reads were cleaned with FastP using the script `scripts/fastp_clean.sh` (executed with `scripts/run_fastp_clean.sh`.)

### Metadata
The metadata was formatted with the scripts `scripts/prepare_desjardins_metadata.xsh` and `scripts/prepare_metadata_ashton.ipynb`. 
It adds the sample ID (SRS accession) to the original metadata tables and standardizes the format. It uses 
`config/Desjardins_Supplemental_Table_S1.csv`, `config/Ashton_Supplemental_Table_S1.csv`, and `config/reads_table.csv` 
files in both `Crypto_Desjardins/` and `Crypto_Ashton/`. It was run with the environment `misc/sra-tools.yaml`. 
It creates the tables `config/metadata.csv` in both `Crypto_Desjardins/` and `Crypto_Ashton/`. 

### Loci table
A table with genes of interest to add to the plots was created with the following information: 
The gene IDs of the centromere adjacent genes were taken from Janbon 2014. 
The MAT loci gene IDs (everything between SXI1 and STE12) were taken from the VNI genome GFF. 
The rRNA genes (with the tags level1 ncRNA and level2 rRNA) from the VNI annotation.    
  * `config/loci.csv`

### RepBase
The RepBase database of consensus repetitive sequences was downloaded with the following steps to create: 
  * `config/RepBase.fasta`  
```
wget https://www.girinst.org/server/RepBase/protected/RepBase29.01.fasta.tar.gz
tar -xvzf RepBase29.01.fasta.tar.gz
cat RepBase29.01.fasta/*.ref > RepBase.fasta
cat RepBase29.01.fasta/appendix/*.ref >> RepBase.fasta
rm -rf RepBase29.01.fasta/ RepBase29.01.fasta.tar.gz 
```

## Workflow execution

1) Run the analysis workflow for each dataset:
```
conda activate snakemake
snakemake --profile Crypto_Desjardins/config/default
snakemake --profile Crypto_Ashton/config/default
``` 
See log files:  
  `Crypto_Desjardins/logs/weavepop_all.log`  
  `Crypto_Ashton/logs/weavepop_all.log`  

2) Run the join datasets workflow:

```
snakemake --profile Crypto_Desjardins_Ashton/config/default
```
See log file:  
`Crypto_Desjardins_Ashton/logs/weavepop_all.log`  

3) Rerun the workflow excluding non-haploid (based on the `explore_depth.qmd` analysis described below)
 and low quality samples. 

See log files:  
  `Crypto_Desjardins/logs/weavepop_exclude.log`  
  `Crypto_Ashton/logs/weavepop_exclude.log`  

4) Rerun the join datasets workflow.  

See log file:  
`Crypto_Desjardins_Ashton/logs/weavepop_exclude.log`  


## Analyses

```
analyses/
├── data # Raw data and data processed by the scripts here
├── misc # `quarto.yaml` environment definition
├── notebooks # Rendered versions of Quarto documents
├── results # Results of the analyses
├── scripts # Scripts, Jupyter notebooks, and Quarto documents
```

The working directory of all the Quarto documents is `analyses/`.  
They are rendered with the command: `quarto render analyses/scripts/<name>.qmd`.  

The script `analyses/scripts/metadata_colors.R` creates color palettes for the metadata to use in the plots.

All the input files used in the analyses come from the input or results of WeavePop in `Crypto_Desjardins`, 
`Crypto_Ashton` or `Crypto_Desjardins_Ashton` or from the following external data.  

External data:  
* `data/raw/CryptoDiversity_Desjardins_Tree.tre`: From [CryptoDiversity_Tree_Info](https://github.com/magwenelab/CryptoDiversity_Tree_Info/blob/main/CryptoDiversity_Desjardins_Tree.tre)
* `data/raw/2017.06.09.all_ours_and_desj.snp_sites.mod.fa.cln.tree`: Received from Philip Ashton (Dec 5 2024).


| Analysis | Script <br /> `scripts/` | Output | Description |
|-----------------|-----------------|-----------------| -----------------|
| Explore Depth Profile of all Samples | `explore_depth.qmd` | `results/tables/ploidy.tsv` | Explore the depth plots to identify putative non-haploid samples to exclude from the analyses. |
| Metadata | `metadata.ipynb` | `data/processed/metadata_all_H99_complete.csv`<br />  `data/processed/metadata_ashton_desj_all_weavepop_final_H99.csv`<br />   `data/processed/metadata_ashton_desj_vni_weavepop_final_H99.csv`  | Create new metadata tables to add the VNI subdivision information from the Ashton study to the Desjardins samples and remove the samples excluded by ploidy or quality. |
| Tree building | `merge_trees.qmd` | `data/processed/tree_ashton.newick`<br />  `data/processed/tree_desjardins.newick`<br />  `data/processed/tree_merged.newick`<br /> Plots in `results/trees/` | Merge the trees of the Ashton and Desjardins datasets. |
| Discover aneuploidies | `aneuploidies.qmd` | `results/tables/chromosome_cnv_categories.tsv`<br /> Plots in `results/figs/`| Categorize chromosomes by coverage of CNVs |
| Plot duplications in tree |`tree_plot_cnvs.qmd`| Plots in `results/trees_dups/`| Plot the merged tree with a heatmap of duplicated chromosomes.|
|Metrics of quality and variants of final dataset | `summary_per_lineage.qmd` | `data/processed/snp_counts_desjardins.csv`<br /> `data/processed/snp_counts_ashton.csv` <br /> `results/tables/per_lineage_summary_stats.tsv`| Create summary table of  mapping stats, number of variants, and CNVs. |
