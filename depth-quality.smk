configfile: "config.yaml"

import pandas as pd

samplefile=(pd.read_csv(config["sample_file"], sep=","))
samples=list(set(samplefile["sample"]))

rule all:
    input:
        expand("genomes-annotations/{sample}/coverage.regions.bed.gz",sample=samples),        
        expand("genomes-annotations/{sample}/coverage.svg",sample=samples),
        expand("genomes-annotations/{sample}/coverage.txt",sample=samples)

rule mosdepth:
    input:
        "genomes-annotations/{sample}/snps.bam"
    output:
        "genomes-annotations/{sample}/coverage.regions.bed.gz"
    params:
        window = config["mosdepth_window"]
    conda: 
        "depth.yaml"
    threads:
       config["mosdepth_threads"]     
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
        "coverage.R"

rule mapq:
    input:
        "genomes-annotations/{sample}/snps.bam"
    output:
        "genomes-annotations/{sample}/mapq.tsv"
    conda: 
        "depth.yaml"   
    log:
        "logs/mapq/{sample}.log"
    shell:
        "samtools stats {input} | grep ^MAPQ | cut -f 2- > {output} "
        "&> {log}"

rule mapq_plot:
    input:
        "genomes-annotations/{sample}/mapq.tsv"
    output:
        "genomes-annotations/{sample}/mapq-count.svg"
    log:
        "logs/mapq-count/{sample}.log"
    script:
        "mapq-count.R"

rule mpileup:
    input:
        "genomes-annotations/{sample}/snps.bam"
    output:
        "genomes-annotations/{sample}/snps.pileup"
    conda: 
        "depth.yaml"   
    log:
        "logs/mapq/{sample}.log"
    shell:
        "samtools mpileup --output-extra MAPQ {input} > {output}"
        "&> {log}"

#samtools stats genomes-annotations/SRS404442/snps.bam > snps.stats
#plot-bamstats -p SRS404442/ snps.stats