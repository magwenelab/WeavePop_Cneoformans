configfile: "config.yaml"

import pandas as pd

samplefile=(pd.read_csv(config["sample_file"], sep=","))
samples=list(set(samplefile["sample"]))

rule all:
    input:
        expand("genomes-annotations/{sample}/coverage.svg",sample=samples),
        expand("genomes-annotations/{sample}/coverage.regions.bed.gz",sample=samples),
        expand("genomes-annotations/{sample}/coverage_good.regions.bed.gz",sample=samples),
        expand("genomes-annotations/{sample}/mapq_distribution.svg",sample=samples),
        expand("genomes-annotations/{sample}/cov_distribution.svg",sample=samples),
        expand("genomes-annotations/{sample}/bamstats", sample=samples),
        "results/mapped_reads.svg",
        expand("genomes-annotations/{sample}/mapq_window.bed", sample=samples),
        expand("genomes-annotations/{sample}/mapq.svg", sample=samples),
        expand("genomes-annotations/{sample}/annotation.gff", sample=samples),
        "results/proprotional_coverage_good.csv",
        "results/cov_good_all.svg",
        "results/cov_prop_median_good.svg",
        "results/cov_prop_mean_good.svg",
        "results/proprotional_coverage_raw.csv",
        "results/cov_raw_all.svg",
        "results/cov_prop_median_raw.svg",
        "results/cov_prop_mean_raw.svg"

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
        "genomes-annotations/{sample}/coverage.regions.bed.gz",
        "genomes-annotations/{sample}/coverage_good.regions.bed.gz",
        "chromosome_names.csv",
        "results/loci_interest.tsv"
    output:
        "genomes-annotations/{sample}/coverage.svg",
        "genomes-annotations/{sample}/coverage_stats.svg",
        "genomes-annotations/{sample}/coverage_raw.csv",
        "genomes-annotations/{sample}/coverage_good.csv",
        "genomes-annotations/{sample}/coverage_global_raw.csv",
        "genomes-annotations/{sample}/coverage_global_good.csv"
    log:
        "logs/coverage/{sample}.log"
    script:
        "scripts/coverage.R"

rule cat_stats:
    input:
        r = expand("genomes-annotations/{sample}/coverage_raw.csv",sample=samples),
        g = expand("genomes-annotations/{sample}/coverage_good.csv",sample=samples),
        gr = expand("genomes-annotations/{sample}/coverage_global_raw.csv",sample=samples),
        gg = expand("genomes-annotations/{sample}/coverage_global_good.csv",sample=samples)
    output:
        allr = "results/coverage_raw.csv",
        allg = "results/coverage_good.csv",
        allgr = "results/coverage_global_raw.csv",
        allgg = "results/coverage_global_good.csv"
    log:
        "logs/coverage/cat_stats.log"
    shell:
        "cat {input.r} > {output.allr} "
        "&& "
        "cat {input.g} > {output.allg} "
        "&& "
        "cat {input.gr} > {output.allgr} "
        "&& "
        "cat {input.gg} > {output.allgg} 2> {log}"

rule coverage_stats_plots:
    input:
        "sample_metadata.csv",
        "results/coverage_global_good.csv",
        "results/coverage_good.csv",
        "results/coverage_global_raw.csv",
        "results/coverage_raw.csv"        
    output:
        "results/proprotional_coverage_good.csv",
        "results/cov_good_all.svg",
        "results/cov_prop_median_good.svg",
        "results/cov_prop_mean_good.svg",
        "results/proprotional_coverage_raw.csv",
        "results/cov_raw_all.svg",
        "results/cov_prop_median_raw.svg",
        "results/cov_prop_mean_raw.svg"
    log:
        "logs/coverage/stats_plot.log"    
    script:
        "scripts/cov_stats_all.R"
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
        "genomes-annotations/{sample}/mapq.csv",
        "chromosome_names.csv"
    output:
        "genomes-annotations/{sample}/mapq_distribution.svg"
    log:
        "logs/mapq-dist/{sample}.log"
    script:
        "scripts/mapq-distribution.R"

rule cov_distribution:
    input:
        "genomes-annotations/{sample}/cov.csv",
        "chromosome_names.csv"
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

rule mapq:
   input:
       "genomes-annotations/{sample}/snps.bam",
       "genomes-annotations/{sample}/coverage.regions.bed.gz"
   output:
       "genomes-annotations/{sample}/mapq.bed",
        "genomes-annotations/{sample}/mapq_window.bed" 
   log:
       "logs/mapq/{sample}.log"
   script:
        "scripts/pileup_mapq.sh"

rule mapq_plot:
    input:
        "genomes-annotations/{sample}/mapq_window.bed",
        "chromosome_names.csv",
        "results/loci_interest.tsv"
    output:
        "genomes-annotations/{sample}/mapq.svg"
    log:
        "logs/mapq_plot/{sample}.log"
    script:
        "scripts/mapq.R"

rule mapqcov2gff:
    input:
        mapqbed = "genomes-annotations/{sample}/mapq_window.bed",
        covbed = "genomes-annotations/{sample}/coverage.regions.bed.gz",
        gff = "genomes-annotations/{sample}/lifted.gff_polished"
    output:
        covmapq = "genomes-annotations/{sample}/mapq_cov_window.bed",
        newgff = "genomes-annotations/{sample}/annotation.gff"
    log: 
        "logs/gff/{sample}.log"
    shell:
        "xonsh scripts/mapqcov2gff.xsh {input.mapqbed} {input.covbed} {input.gff} {output.covmapq} {output.newgff} &> {log}"

