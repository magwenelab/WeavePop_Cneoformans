# Installations
<details>
<summary>Install crypto_div environment </summary>
 
With this `crypto_div.yml` file:
~~~
name: crypto_div
channels:
  - conda-forge
  - defaults
  - bioconda
dependencies:
  - pandas
  - xonsh
  - snakemake
  - snippy
  - liftoff
~~~

Run:
~~~
nohup conda env create -y -f crypto_div.yml &
~~~

When the environment is ready install R:
~~~
conda activate crypto_div
conda install -c r r-essentials
conda deactivate
~~~
And install Graphviz to see Sankemake DAG of jobs in svg
~~~
conda install -c conda-forge graphviz
~~~

</details>

<details>
<summary>Install agat environment </summary>

Run this lines one by one:
~~~
 conda create -n agat
 conda activate agat
 conda install perl-bioperl perl-clone perl-graph perl-lwp-simple perl-carp perl-sort-naturally perl-file-share perl-file-sharedir-install perl-moose perl-yaml perl-lwp-protocol-https -c bioconda
 conda install r-base
 conda install perl-statistics-r -c bioconda
 cpan install bioperl List::MoreUtils Term::ProgressBar
 git clone https://github.com/NBISweden/AGAT.git
 perl Makefile.PL 
 make
 make test
 make install
 conda deactivate
 ~~~
</details>

 # Pipeline modifications
The working directory is `/analysis/czirion/CryptoDiversity`
* `get-lineage-of-samples.xsh`: Adds the SRS codes to the `Desjardins_Supplemental_Table_S1.csv` and puts it in the file 'sample_metadata.csv'.  
* `get-references.xsh`: Uses `lineage_references.csv` and `sample_metadata.csv` to create `sample_reference.csv`, that has each sample name, fastq filenames, lineage, and corresponding reference assembly filename.  
* `fastq2snippy.smk`: Is the Snakefile to run the pipeline, it uses the `config.yaml` file.   For the moment it does the following:  
  * Runs the script `fastq-combiner.xsh` for each sample in `read_pair_table.csv`. This concatenates all `_1.fastq` of one sample into only one file named `<SRS-accession>_1.fq.gz` and compresses it and does the same for `_2.fastq`.  
  * Runs **snippy** for each sample.
* `annotate_references.smk`: Is a Snakefile to lift over annotations from `Reference_Genomes/FungiDB-53_CneoformansH99_PMM.gff` into the four lineages genomes. It currently works with:  
  ` snakemake --snakefile annotate_references.smk --cores 1 --use-conda --conda-frontend conda -p`:  
      ⚠️ `--cores 1` is because there is a problem if liftoff runs in parallel because the different jobs try to create `FungiDB-53_CneoformansH99_PMM.gff_db` at the same time and that is not cool.    
      ⚠️ `--conda-frontend conda` because it cannot use mamba, which is the default.  
      ❔  It makes `{lineage}_liftoff.agat.log` files out of nowhere, they are not specified in the snakefile.  
      ⏰ Pending: Make the name of the `Reference_Genomes` directory configurable in the `config.yaml`  
      ⏰ Pending: Merge into main workflow.  
