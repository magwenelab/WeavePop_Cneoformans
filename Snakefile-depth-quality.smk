configfile: "config.yaml"

import pandas as pd

samplefile=(pd.read_csv(config["sample_file"], sep=","))
samples=list(set(samplefile["sample"]))

rule all:
    input:
        expand("genomes-annotations/{sample}/coverage.svg",sample=samples),
        expand("genomes-annotations/{sample}/coverage.txt",sample=samples),
        expand("genomes-annotations/{sample}/mapq_distribution.svg",sample=samples),
        expand("genomes-annotations/{sample}/cov_distribution.svg",sample=samples),
        expand("genomes-annotations/{sample}/bamstats", sample=samples)

rule mosdepth:
    input:
        "genomes-annotations/{sample}/snps.bam"
    output:
        "genomes-annotations/{sample}/coverage.regions.bed.gz"
    params:
        window = config["mosdepth_window"]
    conda: 
        "envs/depth.yaml"
    threads:
       config["threads_mosdepth"]     
    log:
        "logs/mosdepth/{sample}.log"
    shell:
        "mosdepth -n --by {params.window} -t {threads} "
        "genomes-annotations/{wildcards.sample}/coverage {input} "
        "&> {log}"

rule coverage_plot:
    input:
        "genomes-annotations/{sample}/coverage.regions.bed.gz"
    output:
        "genomes-annotations/{sample}/coverage.txt",
        "genomes-annotations/{sample}/coverage.svg"
    log:
        "logs/coverage/{sample}.log"
    script:
        "scripts/coverage.R"

rule samtools_stats:
    input:
        "genomes-annotations/{sample}/snps.bam",
        "genomes-annotations/{sample}/ref.fa"
    output:
        mapq = "genomes-annotations/{sample}/mapq.tsv",
        cov = "genomes-annotations/{sample}/cov.tsv"
    conda: 
        "envs/depth.yaml"
    log:
        "logs/stats/{sample}.log"
    script:
        "scripts/samtools-stats.xsh {wildcards.sample} &> {log}"

rule mapq_distribution:
    input:
        "genomes-annotations/{sample}/mapq.tsv"
    output:
        "genomes-annotations/{sample}/mapq_distribution.svg"
    log:
        "logs/mapq-dist/{sample}.log"
    script:
        "scripts/mapq-distribution.R"

rule cov_distribution:
    input:
        "genomes-annotations/{sample}/cov.tsv"
    output:
        "genomes-annotations/{sample}/cov_distribution.svg"
    log:
        "logs/cov-dist/{sample}.log"
    script:
        "scripts/coverage-distribution.R"        

#rule mpileup:
#    input:
#        "genomes-annotations/{sample}/snps.bam"
#    output:
#        "genomes-annotations/{sample}/snps.pileup"
#    conda: 
#        "depth.yaml"   
#    log:
#        "logs/mapq/{sample}.log"
#    shell:
#        "samtools mpileup --output-extra MAPQ {input} > {output}"
#        "&> {log}"

rule bamstats:
    input:
        "genomes-annotations/{sample}/snps.bam"
    output:
        "genomes-annotations/{sample}/snps.bam.stats"
    conda:
        "envs/depth.yaml"
    shell:
        "samtools stats {input} > {output} "

rule plot_bamstats:
    input:
        "genomes-annotations/{sample}/snps.bam.stats"
    output:
        directory("genomes-annotations/{sample}/bamstats")
    conda:
        "envs/depth.yaml"
    shell:
        "plot-bamstats -p bamstats {input}"
