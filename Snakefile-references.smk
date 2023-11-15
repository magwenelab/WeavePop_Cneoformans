configfile: "config.yaml"

import pandas as pd

lin_ref_table = (pd.read_csv(config["lineage_reference_file"], sep=","))
LINS=list(lin_ref_table["Group"])
lin_ref_table.set_index('Group', inplace=True)

REFDIR = str(config["reference_directory"])
REF_FASTA = REFDIR + str(config["reference_fasta"])
REF_GFF = REFDIR + str(config["reference_gff"])

rule all:
    input:
        expand(REFDIR + "{lineage}_liftoff.gff_polished",lineage=LINS),
        expand(REFDIR + "{lineage}_predicted_proteins.fa",lineage=LINS),
        expand(REFDIR + "{lineage}_predicted_cds.fa",lineage=LINS),
        expand(REFDIR + "{lineage}_protein_list.txt",lineage=LINS),
        "protein_list.txt"

rule ref2ref_liftoff:
    input:
        target_refs = lambda wildcards: (REFDIR +  lin_ref_table.loc[wildcards.lineage, 'File']),
        fasta = REF_FASTA,
        gff = REF_GFF
    output:
        lin_gff = REFDIR + "{lineage}_liftoff.gff_polished",
        lin_fasta = REFDIR + "{lineage}.fasta"
    threads: config["threads_liftoff"] 
    log:
        "logs/references/{lineage}_ref_liftoff.log"
    params:
        refdir = REFDIR    
    shell:
        "liftoff "
        "-g {input.gff} "
        "-polish "
        "-o {params.refdir}/{wildcards.lineage}_liftoff.gff "
        "-dir {params.refdir}/{wildcards.lineage}_intermediate_files "
        "-u {params.refdir}/{wildcards.lineage}_unmapped_features.txt "
        "-p {threads} "
        "{input.target_refs} {input.fasta} "
        "&> {log} "
        "&& cp {input.target_refs} {output.lin_fasta}"

rule ref2ref_agat:
    input: 
        lin_liftoff = REFDIR + "{lineage}_liftoff.gff_polished",
        lin_fasta = REFDIR + "{lineage}.fasta"
    output:
        cds = REFDIR + "{lineage}_predicted_cds.fa",
        prots = REFDIR + "{lineage}_predicted_proteins.fa"
    conda:
        "envs/agat.yaml"
    log:
        cds = "logs/references/{lineage}_ref_agat_cds.log",   
        prots = "logs/references/{lineage}_ref_agat_prots.log"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.lin_liftoff} "
        "-f {input.lin_fasta} "
        "-o {output.cds} "
        "&> {log.cds} "
        " && "
        "agat_sp_extract_sequences.pl "
        "-g {input.lin_liftoff} "
        "-f {input.lin_fasta} "
        "-o {output.prots} "
        "-p &> {log.prots}"  

rule protein_list:
    input:
        fasta = REFDIR + "{lineage}_predicted_proteins.fa"
    output:
        list = REFDIR + "{lineage}_protein_list.txt"
    shell:
        "seqkit seq -n -i {input.fasta} > {output.list}"

rule cat_lists:
    input: 
        expand(REFDIR + "{lineage}_protein_list.txt", lineage=LINS)
    output:
        "protein_list.txt"
    shell:
        "cat {input} | sort | uniq > {output}"

rule unmapped_features:
    input:
        expand(REFDIR + "{lineage}_unmapped_features.txt", lineage=LINS)   
    output: 
        "unmapped_features.csv"
    shell:
        # cat ../lineage_references.csv | tail -n4 | cut -d',' -f1 | while read lineage
        #   do 
        #    sed "s/$/,\\${lineage}/" ${lineage}_unmapped_features.txt
        #   done > unmapped_features.csv
         
rule unmapped_count:
    R:
    #unmapped<-read.csv("reference_genomes/unmapped_features.csv", header = FALSE, col.names = c("gene", "lineage"))
    #unmapped_count <- unmapped %>% group_by(lineage) %>% count()
    #write.csv(unmapped_count, "reference_genomes/unmapped_count.csv", row.names = FALSE)