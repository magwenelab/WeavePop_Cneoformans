configfile: "config.yaml"

import pandas as pd

samplefile=(pd.read_csv(config["sample_file"], sep=","))
samples=list(set(samplefile["sample"]))

rule all:
    input:
        # expand("genomes-annotations/{sample}/coverage.svg",sample=samples),
        # expand("genomes-annotations/{sample}/coverage.txt",sample=samples),
        # expand("genomes-annotations/{sample}/coverage_good.svg",sample=samples),
        # expand("genomes-annotations/{sample}/coverage_good.txt",sample=samples),
        expand("genomes-annotations/{sample}/coverage.regions.bed.gz",sample=samples),
        expand("genomes-annotations/{sample}/coverage_good.regions.bed.gz",sample=samples),
        expand("genomes-annotations/{sample}/mapq_distribution.svg",sample=samples),
        expand("genomes-annotations/{sample}/cov_distribution.svg",sample=samples),
        expand("genomes-annotations/{sample}/bamstats", sample=samples),
        "results/mapped_reads.svg"

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

rule mosdepth_good:
    input:
        "genomes-annotations/{sample}/snps.bam"
    output:
        "genomes-annotations/{sample}/coverage_good.regions.bed.gz"
    params:
        window = config["mosdepth_window"],
        min_mapq = config["mosdepth_min_mapq"]
    conda: 
        "envs/depth.yaml"
    threads:
       config["threads_mosdepth"]     
    log:
        "logs/mosdepth_good/{sample}.log"
    shell:
        "mosdepth -n --by {params.window} --mapq {params.min_mapq} -t {threads} "
        "genomes-annotations/{wildcards.sample}/coverage_good {input} "
        "&> {log}"

rule coverage_plot:
    input:
        "genomes-annotations/{sample}/coverage.regions.bed.gz"
    output:
        "genomes-annotations/{sample}/coverage.txt",
        "genomes-annotations/{sample}/coverage.svg",
        "genomes-annotations/{sample}/coverage_stats.svg"
    log:
        "logs/coverage/{sample}.log"
    script:
        "scripts/coverage.R"

rule coverage_good_plot:
    input:
        "genomes-annotations/{sample}/coverage_good.regions.bed.gz"
    output:
        "genomes-annotations/{sample}/coverage_good.txt",
        "genomes-annotations/{sample}/coverage_good.svg",
        "genomes-annotations/{sample}/coverage_good_stats.svg"
    log:
        "logs/coverage_good/{sample}.log"
    script:
        "scripts/coverage.R"

rule samtools_stats:
    input:
        bam = "genomes-annotations/{sample}/snps.bam",
        ref = "genomes-annotations/{sample}/ref.fa"
    output:
        mapq = "genomes-annotations/{sample}/mapq.csv",
        cov = "genomes-annotations/{sample}/cov.csv"
    log:
        "logs/stats/{sample}.log"
    shell:
        "xonsh scripts/samtools-stats.xsh {wildcards.sample} {input.bam} {input.ref} {output.mapq} {output.cov} &> {log}"

rule mapq_distribution:
    input:
        "genomes-annotations/{sample}/mapq.csv"
    output:
        "genomes-annotations/{sample}/mapq_distribution.svg"
    log:
        "logs/mapq-dist/{sample}.log"
    script:
        "scripts/mapq-distribution.R"

rule cov_distribution:
    input:
        "genomes-annotations/{sample}/cov.csv"
    output:
        "genomes-annotations/{sample}/cov_distribution.svg"
    log:
        "logs/cov-dist/{sample}.log"
    script:
        "scripts/coverage-distribution.R"        
rule bamstats:
    input:
        "genomes-annotations/{sample}/snps.bam"
    output:
        "genomes-annotations/{sample}/snps.bam.stats"
    conda:
        "envs/depth.yaml"
    log:
        "logs/bamstats/{sample}.log"
    shell:
        "samtools stats {input} 1> {output} 2> {log}"

rule plot_bamstats:
    input:
        "genomes-annotations/{sample}/snps.bam.stats"
    output:
        directory("genomes-annotations/{sample}/bamstats")
    conda:
        "envs/depth.yaml"
    log:
        "logs/plot-bamstats/{sample}.log"
    shell:
        "plot-bamstats -p {output}/ {input} &> {log}"

rule unmapped_edit:
    input:
        "genomes-annotations/{sample}/snps.bam.stats" 
    output: 
        temp("genomes-annotations/{sample}/mapping_stats.txt")
    shell:
        "grep reads {input} | cut -d'#' -f1 | cut -f 2- | grep . > {output} "
        " && "
        'sed -i "s/$/:\\{wildcards.sample}/" {output}'

rule unmapped:
    input:
        expand("genomes-annotations/{sample}/mapping_stats.txt", sample=samples)   
    output: 
        "results/mapping_stats.txt"
    shell:
       'cat {input} > {output}'  

rule unmapped_plot:
    input:
        "results/mapping_stats.txt",
        "sample_metadata.csv"
    output:
        "results/mapped_reads.svg"
    log:
        "logs/stats/mapped.log"
    script:
        "scripts/mapped_reads.R"

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