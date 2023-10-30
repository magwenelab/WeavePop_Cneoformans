configfile: "config.yaml"

import pandas as pd
from glob import glob

samplefile=(pd.read_csv(config["sample_file"], sep=","))
samples=list(set(samplefile["sample"]))
ref_table = (pd.read_csv(config["sample_reference_file"], sep=","))
ref_table.set_index('sample', inplace=True)
REFDIR = str(config["reference_directory"])

protlist=(pd.read_csv("protein_list.txt", sep=",", header = None, names = ['protein']))
proteins=list(protlist["protein"])

#INPUT_FILE = 'proteins/{sample}_{protein}.fa'
#INP = glob_wildcards(INPUT_FILE).sample

rule all:
    input:
        expand("genomes-annotations/{sample}/snps.consensus.fa",sample=samples),
        expand("genomes-annotations/{sample}/lifted.gff_polished", sample=samples),
        expand("genomes-annotations/{sample}/predicted_cds.fa",sample=samples),
        expand("genomes-annotations/{sample}/predicted_proteins.fa",sample=samples),
        expand("cds/{sample}.done",sample=samples),
        expand("proteins/{sample}.done",sample=samples),
        #expand("cds/{protein}.fa", protein=proteins),
        #expand("proteins/{protein}.fa", protein=proteins)

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
        "genomes-annotations/{sample}/snps.consensus.fa"
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

rule by_cds:
    input:
        fasta = "genomes-annotations/{sample}/predicted_cds.fa",
        list = "protein_list.txt",
        idx = "genomes-annotations/{sample}/predicted_cds.fa.fai"
    output:
        temp("cds/{sample}.done")
    conda:
        "agat.yaml"
    log:
        "logs/cds/{sample}.log"    
    shell:
        "cat {input.list} | "
        "while read line;  do "
        "seqkit faidx {input.fasta} $line | "
        "seqkit replace -p '($)' -r ' sample={wildcards.sample}' > cds/{wildcards.sample}_$line.fa; done &> {log}"
        "&& touch {output}" 

rule by_protein:
    input:
        fasta = "genomes-annotations/{sample}/predicted_proteins.fa",
        list = "protein_list.txt",
        idx = "genomes-annotations/{sample}/predicted_proteins.fa.fai"
    output:
        temp("proteins/{sample}.done")
    conda:
        "agat.yaml"
    log:
        "logs/proteins/{sample}.log"   
    shell:
        "cat {input.list} | "
        "while read line; do "
        "seqkit faidx {input.fasta} $line | "
        "seqkit replace -p '($)' -r ' sample={wildcards.sample}' > proteins/{wildcards.sample}_$line.fa; done &> {log}"
        "&& touch {output}"

#rule concatenate_prots:
#    input:
#        expand(INPUT_FILE, sample=INP, protein="{protein}"),
#    output:
#        'proteins/{protein}.fa'
#    shell:
#        "cat {input} > {output} "
#rule concatenate_cds:
#    input:
#        expand(INPUT_FILE, sample=INP, protein="{protein}"),
#    output:
#        'cds/{protein}.fa'
#    shell:
#        "cat {input} > {output} "
