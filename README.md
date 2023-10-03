# Crypto_Diversity_Pipeline

Pipeline for genome assembly and analysis of data from Desjardins et al. 2017


## Requirements

* Python and Conda (miniconda is my preference)
* Xonsh -- https://xon.sh/
* Snakemake -- https://snakemake.github.io/
* NCBI Entrez Utilities (E-utilities) command line tools -- https://www.ncbi.nlm.nih.gov/books/NBK25501/
* NCBI SRA Tools -- https://github.com/ncbi/sra-tools


The following are most easily installed via Conda. AGAT in particular seems to clash with other tools so I install it in its own environment.

* Python modules -- Pandas
* R with tidyverse meta-package
* Snippy -- https://github.com/tseemann/snippy
* Liftoff -- https://github.com/agshumate/Liftoff
* AGAT -- https://github.com/NBISweden/AGAT
  


## Overview

Scripts to be run in this order:

1. `CryptoDiversity-Retrieve.xsh` -- given an NCBI BioProject ID, identifies all the BioSamples associated with that project and downloads each into a folder called `Samples/${SRSID}` where `SRSID` are SRA SRS numbers
  
    * Files produced:
      
      * "{project_id}_samples.txt" -- a tab-delimited text file listing all the samples associated with the BioProject.  Three columns: SAM, Name, SRS  where SAM = BioSample Sample ID , Name = common name, SRS = SRA sequence read identifier
      
      * "{project_id}_SRStoSRR.txt" -- tab-delimited text file giving mapping between SRS (Samples) IDs and SRR (Runs) IDs. A given SRS can have multiple associated SRRs. 
      
      * Each `Samples/${SRSID}` directory will contain a number of SRR subdirectories of the form `Samples/${SRSID}/{SRRID}` where SRA runs are "prefetched" using the SRA Tools `prefetch` command.


2. `fastq-unpacker.xsh` -- turn each of the "prefetched" run files (produced in the prior step) into FASTQ files, written to `./fastqs` directory

3. `build-read-pair-tables.xsh` -- create CSV files for paired and unpaired fastqs
   * Files produced:
  
     * `read_pair_table.csv` -- Columns are SRS ID, SRR ID, read pair 1, read pair 2, total size in bytes of read pair

     * `largest_read_pair_table.csv` -- A filtered version of `read_pair_table.csv`, giving the single largest pair per sample
  
     * `unpaired_fastqs.csv` -- Columns are SRS ID, SRR ID, fastq file(s)


4. `snippy-builder.xsh` -- runs Snippy pipeline on each of the read pairs in `largest_read_pair_table.csv`.  This step is a key one for optimization (perhaps moving to Snakemake so can run in parallel) and improvement (e.g. don't just use largest read pair, but all read pairs)
    
    * Generates reference alignments and SNP calls for each sample

    * Improvements to make: 
      * Don't just use largest read pair, but all read pairs
      * Use appropriate reference for each strain, not the same reference (see supplemental tables from Desjardins et al. for lineage date).
      * Maybe move this to Snakemake?


5. `Snakefile.liftoff` -- should be run from w/in the `Samples/` directory (e.g. create soft link).  Runs annotation liftover via liftoff, generates predicted CDS and protein sequences for each gene via AGAT, for each gene creates files with all predicted CDS and AA seqs for that gene (across the 384 samples)





## Notes

Here's a list of potentially suitable reference strains for each of the C. neoformans lineages

* VNI -- H99, GCA_011801205.1
* VNII -- VNII strain collected from cockatoo excrement , GCA_022832995.1
* VNBI -- Ftc146-1, GCA_010065285.1
* VNBII -- Bt81, GCA_023650555.1

Downside of using the above Genbank assemblies, as opposed to using FungiDB sources, is that not all of the Genbank assemblies have included annotation or where they do have annotation, not using the same scheme.  This would complicate the step of pulling all the sequences of a given ORF into one file.

Alternative might be to use liftoff to "lift over" annotation from reference H99 genome to each of the lineage specific references prior to analysis.

