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
        "protein_list.txt",
        REFDIR + "references_unmapped_features.csv",
        REFDIR + "references_unmapped_count.csv"

rule features:
    input:
        REF_GFF
    output:
        REFDIR + "features.txt"
    shell:
        "grep -v '#' {input} | cut -f 3 | sort | uniq > {output}"

rule ref2ref_liftoff:
    input:
        target_refs = REFDIR + "{lineage}.fasta",
        fasta = REF_FASTA,
        gff = REF_GFF,
        features = REFDIR + "features.txt"
    output:
        lin_gff = REFDIR + "{lineage}_liftoff.gff_polished",
    threads: config["threads_liftoff"] 
    log:
        "logs/references/{lineage}_ref_liftoff.log"
    params:
        refdir = REFDIR    
    shell:
        "liftoff "
        "-g {input.gff} "
        "-polish "
        "-f {input.features} "
        "-o {params.refdir}/{wildcards.lineage}_liftoff.gff "
        "-dir {params.refdir}/{wildcards.lineage}_intermediate_files "
        "-u {params.refdir}/{wildcards.lineage}_unmapped_features.txt "
        "-p {threads} "
        "{input.target_refs} {input.fasta} "
        "&> {log} "

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

rule unmapped_features_edit:
    input:
        REFDIR + "{lineage}_unmapped_features.txt"   
    output: 
        temp(REFDIR + "{lineage}_unmapped_features.csv")
    shell:
       'sed "s/$/,\\{wildcards.lineage}/" {input} > {output}'

rule unmapped_features:
    input:
        expand(REFDIR + "{lineage}_unmapped_features.csv", lineage=LINS)   
    output: 
        REFDIR + "references_unmapped_features.csv"
    shell:
       'cat {input} > {output}'         

rule unmapped_count:
    input:
        REFDIR + "references_unmapped_features.csv"
    output:
        REFDIR + "references_unmapped_count.csv"
    script:
        "scripts/count_reference_unmapped.R"