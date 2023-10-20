configfile: "config.yaml"

import pandas as pd
samplefile=(pd.read_csv(config["sample_file"], sep=","))
samples=list(set(samplefile["sample"]))
ref_table = (pd.read_csv(config["sample_reference_file"], sep=","))
ref_table.set_index('sample', inplace=True)

rule all:
    input:
        expand("fastq_combined/{sample}_1.fq.gz",sample=samples),
        expand("fastq_combined/{sample}_2.fq.gz",sample=samples)

rule combine_fastq:
    input:
        config["sample_file"]
    output:
        "fastq_combined/{sample}_1.fq.gz",
        "fastq_combined/{sample}_2.fq.gz"
    shell:
        "xonsh fastq-combiner.xsh {wildcards.sample} {input} fastqs/ fastq_combined/"

rule snippy:
    input:
        "fastq_combined/{sample}_1.fq.gz",
        "fastq_combined/{sample}_2.fq.gz",
        #refdir = config["reference_directory"]
    params:
        ref = lambda wildcards: ref_table.loc[wildcards.sample, 'refgenome'],
        file1 = lambda wildcards: ref_table.loc[wildcards.sample, 'file1'],
        file2 = lambda wildcards: ref_table.loc[wildcards.sample, 'file2'] 
    output:
        directory("snippy-analysis/{sample}")
    threads: 12   
    log:
        "logs/snippy/{sample}.log" 
    shell:
        "snippy --outdir {output} --cpus {threads} "
        "--ref Reference_Genomes/{params.ref} "
        "--R1 fastq_combined/{params.file1} --R2 fastq_combined/{params.file2} &> {log}"
