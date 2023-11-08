configfile: "config.yaml"

import pandas as pd

samplefile=(pd.read_csv(config["sample_file"], sep=","))
samples=list(set(samplefile["sample"]))
ref_table = (pd.read_csv(config["sample_reference_file"], sep=","))
ref_table.set_index('sample', inplace=True)
REFDIR = str(config["reference_directory"])

protlist=(pd.read_csv("protein_list.txt", sep=",", header = None, names = ['protein']))
proteins=list(protlist["protein"])

rule all:
    input:
        expand("genomes-annotations/{sample}/snps.consensus.fa",sample=samples),
        expand("genomes-annotations/{sample}/lifted.gff_polished", sample=samples),
        expand("genomes-annotations/{sample}/predicted_cds.fa",sample=samples),
        expand("genomes-annotations/{sample}/predicted_proteins.fa",sample=samples),
        expand("by_cds/{protein}.fa", protein=proteins),
        expand("by_protein/{protein}.fa", protein=proteins),   
        expand("genomes-annotations/{sample}/coverage.svg",sample=samples),
        expand("genomes-annotations/{sample}/coverage.txt",sample=samples),
        expand("genomes-annotations/{sample}/mapq.tsv",sample=samples),
        expand("genomes-annotations/{sample}/mapq-count.svg",sample=samples)

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
        "xonsh fastq-combiner.xsh {wildcards.sample} {input.readtab} {input.fqdir} fastq_combined/  &> {log}"

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
        "genomes-annotations/{sample}/snps.consensus.fa"
    params:
        refgff = lambda wildcards:(REFDIR + ref_table.loc[wildcards.sample, 'lineage'] + "_liftoff.gff_polished"),
        refgenome = lambda wildcards:(REFDIR + ref_table.loc[wildcards.sample, 'refgenome'])
    output:
        "genomes-annotations/{sample}/lifted.gff",        
        "genomes-annotations/{sample}/lifted.gff_polished"
    threads: config["threads_liftoff"]
    log:
        "logs/liftoff/{sample}.log" 
    shell:
        "liftoff "
        "-g {params.refgff} "
        "-polish "
        "-dir genomes-annotations/{wildcards.sample}/intermediate_files "
        "-u genomes-annotations/{wildcards.sample}/unmapped_features.txt "
        "-o genomes-annotations/{wildcards.sample}/lifted.gff "
        "-p {threads} "
        "{input} "
        "{params.refgenome} &> {log}"

rule agat:
    input:
        "genomes-annotations/{sample}/lifted.gff_polished"
    output:
        cds = "genomes-annotations/{sample}/predicted_cds.fa",
        prots = "genomes-annotations/{sample}/predicted_proteins.fa"
    conda:
        "agat.yaml"
    log: 
        cds = "logs/agat/{sample}_cds.log",
        prots = "logs/agat/{sample}_prots.log"
    shell:
        "seqkit seq -w 70 genomes-annotations/{wildcards.sample}/snps.consensus.fa > "
        "genomes-annotations/{wildcards.sample}/snps.consensus.wrapped.fa && "
        "agat_sp_extract_sequences.pl "
        "-g {input} " 
        "-f genomes-annotations/{wildcards.sample}/snps.consensus.wrapped.fa "
        "-o {output.cds} "
        "&> {log.cds} "
        " && "
        "agat_sp_extract_sequences.pl "
        "-g {input} " 
        "-f genomes-annotations/{wildcards.sample}/snps.consensus.wrapped.fa "
        "-o {output.prots} "
        "-p  &> {log.prots}" 

rule index_proteins:
    input:
        "genomes-annotations/{sample}/predicted_proteins.fa"
    output:
        "genomes-annotations/{sample}/predicted_proteins.fa.fai"
    conda:
        "agat.yaml"        
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
        "agat.yaml"        
    log:
        "logs/faidx/{sample}_cds.log"      
    shell:
        "seqkit faidx {input} &> {log}"

rule get_cds:
    input:
        fasta = "genomes-annotations/{sample}/predicted_cds.fa",
        list = "protein_list.txt",
        idx = "genomes-annotations/{sample}/predicted_cds.fa.fai"
    output:
        done = temp("cds/{sample}.done"),
        fas = expand("cds/{{sample}}_{protein}.fa", protein=proteins)
    conda:
        "agat.yaml"
    log:
        "logs/cds/{sample}.log"    
    shell:
        "cat {input.list} | "
        "while read line;  do "
        "seqkit faidx {input.fasta} $line | "
        "seqkit replace -p '($)' -r ' sample={wildcards.sample}' > cds/{wildcards.sample}_$line.fa; done &> {log}"
        "&& touch {output.done}" 

rule get_protein:
    input:
        fasta = "genomes-annotations/{sample}/predicted_proteins.fa",
        list = "protein_list.txt",
        idx = "genomes-annotations/{sample}/predicted_proteins.fa.fai"
    output:
        done = temp("proteins/{sample}.done"),
        fas = expand("proteins/{{sample}}_{protein}.fa", protein=proteins)
    conda:
        "agat.yaml"
    log:
        "logs/proteins/{sample}.log"   
    shell:
        "cat {input.list} | "
        "while read line; do "
        "seqkit faidx {input.fasta} $line | "
        "seqkit replace -p '($)' -r ' sample={wildcards.sample}' > proteins/{wildcards.sample}_$line.fa; done &> {log}"
        "&& touch {output.done}"

rule cat_proteins:
    input:
        fastas = expand("proteins/{sample}_{{protein}}.fa", sample=samples),
        done = expand("proteins/{sample}.done", sample=samples)
    output:
        "by_protein/{protein}.fa"
    shell:
        "cat {input.fastas} > {output}"

rule cat_cds:
    input:
        fastas = expand("cds/{sample}_{{protein}}.fa", sample=samples),
        done = expand("cds/{sample}.done", sample=samples)
    output:
        "by_cds/{protein}.fa"
    shell:
        "cat {input.fastas} > {output}"


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
