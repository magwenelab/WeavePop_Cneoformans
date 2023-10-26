configfile: "config.yaml"

import pandas as pd

samplefile=(pd.read_csv(config["sample_file"], sep=","))
samples=list(set(samplefile["sample"]))
ref_table = (pd.read_csv(config["sample_reference_file"], sep=","))
ref_table.set_index('sample', inplace=True)

rule all:
    input:
        expand("snippy-analysis/{sample}",sample=samples),
        expand("snippy-analysis/{sample}/liftoff/lifted.gff_polished", sample=samples)
rule combine_fastq:
    input:
        readtab = config["sample_file"],
        fqdir = str(config["fastq_dir"])
    output:
        "fastq_combined/{sample}_1.fq.gz",
        "fastq_combined/{sample}_2.fq.gz"
    log:
        "logs/combine_fastq/{sample}.log"  
    shell:
        "xonsh fastq-combiner.xsh {wildcards.sample} {input.readtab} {input.fqdir} fastq_combined/  &> {log}"

rule snippy:
    input:
        "fastq_combined/{sample}_1.fq.gz",
        "fastq_combined/{sample}_2.fq.gz",
        refdir = str(config["reference_directory"])
    params:
        ref = lambda wildcards: ref_table.loc[wildcards.sample, 'refgenome'],
        file1 = lambda wildcards: ref_table.loc[wildcards.sample, 'file1'],
        file2 = lambda wildcards: ref_table.loc[wildcards.sample, 'file2'] 
    output:
        "snippy-analysis/{sample}/snps.consensus.fa"
    threads: config["threads_snippy"]
    log:
        "logs/snippy/{sample}.log" 
    shell:
        "snippy --outdir snippy-analysis/{wildcards.sample} --cpus {threads} "
        "--ref {input.refdir}/{params.ref} "
        "--R1 fastq_combined/{params.file1} --R2 fastq_combined/{params.file2}  --force &> {log}"

rule liftoff:
    input:
        "snippy-analysis/{sample}/snps.consensus.fa"
    params:
        refgff = lambda wildcards:("Reference_Genomes/" + ref_table.loc[wildcards.sample, 'lineage'] + "_liftoff.gff_polished"),
        refgenome = lambda wildcards:("Reference_Genomes/" + ref_table.loc[wildcards.sample, 'refgenome']),
    output:
        "snippy-analysis/{sample}/liftoff/lifted.gff",        
        "snippy-analysis/{sample}/liftoff/lifted.gff_polished"
    log:
        "logs/liftoff/{sample}.log"    
    shell:
        "liftoff "
        "-g {params.refgff} "
        "-polish "
        "-dir snippy-analysis/{wildcards.sample}/liftoff/intermediate_files "
        "-u snippy-analysis/{wildcards.sample}/liftoff/unmapped_features.txt "
        "-o snippy-analysis/{wildcards.sample}/liftoff/lifted.gff "
        "{input} "
        "{params.refgenome} &> {log}"
