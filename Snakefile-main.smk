configfile: "config.yaml"

#subworkflow otherworkflow:
#    workdir:
#        "."
#    snakefile:
#        "./annotate_references.smk"
#    configfile:
#        "./config.yaml"

import pandas as pd

samplefile=(pd.read_csv(config["sample_file"], sep=","))
samples=list(set(samplefile["sample"]))
ref_table = (pd.read_csv(config["sample_reference_file"], sep=","))
ref_table.set_index('sample', inplace=True)
REFDIR = str(config["reference_directory"])

protlist=(pd.read_csv("results/protein_list.txt", sep=",", header = None, names = ['protein']))
proteins=list(protlist["protein"])


rule all:
    input:
        expand("genomes-annotations/{sample}/snps.consensus.fa",sample=samples),
        expand("genomes-annotations/{sample}/lifted.gff_polished", sample=samples),
        expand("genomes-annotations/{sample}/predicted_cds.fa",sample=samples),
        expand("genomes-annotations/{sample}/predicted_proteins.fa",sample=samples),
        #expand("by_cds/{protein}.fa", protein=proteins),
        #expand("by_protein/{protein}.fa", protein=proteins),
        "results/samples_unmapped.png"
rule combine_fastq:
    input:
        readtab = config["sample_file"],
        fqdir = str(config["fastq_dir"])
    output:
        "fastq_combined/{sample}_1.fq.gz",
        "fastq_combined/{sample}_2.fq.gz"
    log:
        "logs/fastq-combiner/{sample}.log"  
    shell:
        "xonsh scripts/fastq-combiner.xsh {wildcards.sample} {input.readtab} {input.fqdir} fastq_combined/  &> {log}"

rule snippy:
    input:
        "fastq_combined/{sample}_1.fq.gz",
        "fastq_combined/{sample}_2.fq.gz"
    params:
        ref = lambda wildcards: (REFDIR + ref_table.loc[wildcards.sample, 'refgenome']),
        file1 = lambda wildcards: ref_table.loc[wildcards.sample, 'file1'],
        file2 = lambda wildcards: ref_table.loc[wildcards.sample, 'file2'] 
    output:
        "genomes-annotations/{sample}/snps.consensus.fa",
        "genomes-annotations/{sample}/snps.bam"
    threads: config["threads_snippy"]
    log:
        "logs/snippy/{sample}.log" 
    shell:
        "snippy --outdir genomes-annotations/{wildcards.sample} "
        "--cpus {threads} "
        "--ref {params.ref} "
        "--R1 fastq_combined/{params.file1} "
        "--R2 fastq_combined/{params.file2} "
        "--force &> {log}"

rule liftoff:
    input:
        target = "genomes-annotations/{sample}/snps.consensus.fa",
        features = REFDIR + "features.txt"
    params:
        refgff = lambda wildcards:(REFDIR + ref_table.loc[wildcards.sample, 'lineage'] + "_liftoff.gff_polished"),
        refgenome = lambda wildcards:(REFDIR + ref_table.loc[wildcards.sample, 'refgenome'])
    output:
        "genomes-annotations/{sample}/lifted.gff",        
        "genomes-annotations/{sample}/lifted.gff_polished",
        "genomes-annotations/{sample}/unmapped_features.txt"
    threads: config["threads_liftoff"]
    log:
        "logs/liftoff/{sample}.log" 
    shell:
        "liftoff "
        "-g {params.refgff} "
        "-polish "
        "-f {input.features} "
        "-dir genomes-annotations/{wildcards.sample}/intermediate_files "
        "-u genomes-annotations/{wildcards.sample}/unmapped_features.txt "
        "-o genomes-annotations/{wildcards.sample}/lifted.gff "
        "-p {threads} "
        "{input.target} "
        "{params.refgenome} &> {log}"

rule agat:
    input:
        gff = "genomes-annotations/{sample}/lifted.gff_polished",
        fa = "genomes-annotations/{sample}/snps.consensus.fa"
    output:
        cds = "genomes-annotations/{sample}/predicted_cds.fa",
        prots = "genomes-annotations/{sample}/predicted_proteins.fa"
    conda:
        "envs/agat.yaml"
    log: 
        cds = "logs/agat/{sample}_cds.log",
        prots = "logs/agat/{sample}_prots.log"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} " 
        "-f {input.fa} "
        "-o {output.cds} "
        "&> {log.cds} "
        " && "
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} " 
        "-f {input.fa} "
        "-o {output.prots} "
        "-p  &> {log.prots}" 

rule index_proteins:
    input:
        "genomes-annotations/{sample}/predicted_proteins.fa"
    output:
        "genomes-annotations/{sample}/predicted_proteins.fa.fai"
    conda:
        "envs/agat.yaml"        
    log:
        "logs/faidx/{sample}_proteins.log"    
    shell:
        "seqkit faidx {input} &> {log}"

rule index_cds:
    input:
        "genomes-annotations/{sample}/predicted_cds.fa"
    output:
        "genomes-annotations/{sample}/predicted_cds.fa.fai"
    conda:
        "envs/agat.yaml"        
    log:
        "logs/faidx/{sample}_cds.log"      
    shell:
        "seqkit faidx {input} &> {log}"

rule get_cds:
    input:
        fasta = "genomes-annotations/{sample}/predicted_cds.fa",
        list = "results/protein_list.txt",
        idx = "genomes-annotations/{sample}/predicted_cds.fa.fai"
    output:
        done = temp("all_cds/{sample}.done"),
        fas = expand("all_cds/{{sample}}_{protein}.fa", protein=proteins)
    conda:
        "envs/agat.yaml"
    log:
        "logs/all_cds/{sample}.log"    
    shell:
        "cat {input.list} | "
        "while read line;  do "
        "seqkit faidx {input.fasta} $line | "
        "seqkit replace -p '($)' -r ' sample={wildcards.sample}' > all_cds/{wildcards.sample}_$line.fa; done &> {log}"
        "&& touch {output.done}" 

rule get_protein:
    input:
        fasta = "genomes-annotations/{sample}/predicted_proteins.fa",
        list = "results/protein_list.txt",
        idx = "genomes-annotations/{sample}/predicted_proteins.fa.fai"
    output:
        done = temp("all_proteins/{sample}.done"),
        fas = expand("all_proteins/{{sample}}_{protein}.fa", protein=proteins)
    conda:
        "envs/agat.yaml"
    log:
        "logs/all_proteins/{sample}.log"   
    shell:
        "cat {input.list} | "
        "while read line; do "
        "seqkit faidx {input.fasta} $line | "
        "seqkit replace -p '($)' -r ' sample={wildcards.sample}' > all_proteins/{wildcards.sample}_$line.fa; done &> {log}"
        "&& touch {output.done}"

rule cat_proteins:
    input:
        fastas = expand("all_proteins/{sample}_{{protein}}.fa", sample=samples),
        done = expand("all_proteins/{sample}.done", sample=samples)
    output:
        "by_protein/{protein}.fa"
    log:
        "logs/bash/{protein}_cat_proteins.log"
    shell:
        "cat {input.fastas} > {output} 2> {log}"

rule cat_cds:
    input:
        fastas = expand("all_cds/{sample}_{{protein}}.fa", sample=samples),
        done = expand("all_cds/{sample}.done", sample=samples)
    output:
        "by_cds/{protein}.fa"
    log:
        "logs/bash/{protein}_cat_cds.log"
    shell:
        "cat {input.fastas} > {output} 2> {log}"

rule unmapped_features_edit:
    input:
        "genomes-annotations/{sample}/unmapped_features.txt"   
    output: 
        temp("genomes-annotations/{sample}/unmapped_features.csv")
    log:
        "logs/bash/{sample}_unmapped_features_edit.log"
    shell:
       'sed "s/$/,\\{wildcards.sample}/" {input} > {output} 2> {log}'

rule unmapped_features:
    input:
        expand("genomes-annotations/{sample}/unmapped_features.csv", sample=samples)   
    output: 
        "results/samples_unmapped_features.csv"
    log:
        "logs/bash/unmapped_features.log"
    shell:
       'cat {input} > {output} 2> {log}'         

rule unmapped_count:
    input:
        "results/samples_unmapped_features.csv",
        REFDIR + "reference_genes.tsv",
        config["sample_reference_file"]
    output:
        "results/samples_unmapped_count.csv",
        "results/samples_unmapped.png"
    log:
        "logs/bash/unmapped_count.log"
    script:
        "scripts/count_sample_unmapped.R"