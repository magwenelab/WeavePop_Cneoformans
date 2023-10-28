configfile: "config.yaml"

import pandas as pd

samplefile=(pd.read_csv(config["sample_file"], sep=","))
samples=list(set(samplefile["sample"]))
ref_table = (pd.read_csv(config["sample_reference_file"], sep=","))
ref_table.set_index('sample', inplace=True)
REFDIR = str(config["reference_directory"])

rule all:
    input:
        expand("genomes-annotations/{sample}/snps.consensus.fa",sample=samples),
        expand("genomes-annotations/{sample}/lifted.gff_polished", sample=samples),
        expand("genomes-annotations/{sample}/predicted_cds.fa",sample=samples),
        expand("genomes-annotations/{sample}/predicted_proteins.fa",sample=samples),
        expand("genomes-annotations/{sample}/predicted_proteins.fa.fai",sample=samples),
        expand("cds/{sample}.done",sample=samples),
        expand("proteins/{sample}.done",sample=samples)

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
        "fastq_combined/{sample}_2.fq.gz",
        refdir = REFDIR
    params:
        ref = lambda wildcards: ref_table.loc[wildcards.sample, 'refgenome'],
        file1 = lambda wildcards: ref_table.loc[wildcards.sample, 'file1'],
        file2 = lambda wildcards: ref_table.loc[wildcards.sample, 'file2'] 
    output:
        "genomes-annotations/{sample}/snps.consensus.fa"
    threads: config["threads_snippy"]
    log:
        "logs/snippy/{sample}.log" 
    shell:
        "snippy --outdir genomes-annotations/{wildcards.sample} "
        "--cpus {threads} "
        "--ref {input.refdir}/{params.ref} "
        "--R1 fastq_combined/{params.file1} "
        "--R2 fastq_combined/{params.file2} "
        "--force &> {log}"

rule liftoff:
    input:
        "genomes-annotations/{sample}/snps.consensus.fa"
    params:
        refgff = lambda wildcards:(REFDIR + ref_table.loc[wildcards.sample, 'lineage'] + "_liftoff.gff_polished"),
        refgenome = lambda wildcards:(REFDIR + ref_table.loc[wildcards.sample, 'refgenome']),
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
    shell:
        "seqkit faidx {input}"

rule sample_list:
    output:
        "samples.txt"
    shell:
        "ls -d genomes-annotations/SRS* > {output}"     

rule by_cds:
    input:
        fasta = "genomes-annotations/{sample}/predicted_cds.fa",
        list = "protein_list.txt"
    output:
        "cds/{sample}.done"
    shell:
        "cat protein_list.txt | "
        "while read line; do "
        "seqkit faidx genomes-annotations/{wildcards.sample}/predicted_cds.fa $line | "
        "seqkit replace -p '($)' -r ' sample={wildcards.sample}' >> cds/$line.fa; done "
        "&& touch {output}" 

rule by_protein:
    input:
        fasta = "genomes-annotations/{sample}/predicted_proteins.fa",
        list = "protein_list.txt"
    output:
        "proteins/{sample}.done"  
    shell:
        "cat protein_list.txt | "
        "while read line; do "
        "seqkit faidx genomes-annotations/{wildcards.sample}/predicted_proteins.fa $line | "
        "seqkit replace -p '($)' -r ' sample={wildcards.sample}' >> proteins/$line.fa; done "
        "&& touch {output}"
