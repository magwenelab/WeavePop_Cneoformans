# Installations
Installing `crypto_div` environment with this `crypto_div.yml` file:
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
  - sra-tools
  - entrez-direct
  - r
  - snippy
  - liftoff
~~~

And running:
~~~
nohup conda env create -y -f crypto_div.yml &
~~~