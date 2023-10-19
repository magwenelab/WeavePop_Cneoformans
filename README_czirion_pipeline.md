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

* 'fastq-combiner.xsh' concatenates all `_1.fastq` of one sample into only one file named `<SRS-accession>_1.fq.gz` and compresses it and does the same for `_2.fastq`.  
* `fastq2snippy.smk` is the Snakefile to run the pipeline, for the moment it runs the previous script for each sample in `read_pair_table.csv`. Using the `config.yaml` file.  
* 'get-lineage-of-samples.xsh' adds the SRS codes to the Desjardines supplemental table and puts it in the file 'sample_metadata.csv'.
* `snippy-builder.xsh` uses the `sample_metadata.csv` and `lineage_references.csv` to run **snippy** in each sample using the appropriate reference genome.  
