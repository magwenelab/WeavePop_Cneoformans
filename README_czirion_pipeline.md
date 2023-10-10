# Installations
+ Install `crypto_div` environment  
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

+ Install `agat` environment  
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

 # Pipeline for snippy-builder
 Bring data to my user:
 ~~~
 cd /analysis/czirion/CryptoDiversity
 ln -s /data/sequence-data/CryptoDiversity/fastqs/ .
 cp /data/sequence-data/CryptoDiversity/unpaired_fastqs.csv .
 cp /data/sequence-data/CryptoDiversity/read_pair_table.csv .
 cp /data/sequence-data/CryptoDiversity/largest_read_pair_table.
 cp /data/sequence-data/CryptoDiversity/snippy-builder.xsh .
~~~
