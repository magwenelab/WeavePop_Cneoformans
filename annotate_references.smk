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
        expand(REFDIR + "{lineage}_predicted_cds.fa",lineage=LINS)

rule ref2ref_liftoff:
    input:
        target_refs = lambda wildcards: (REFDIR +  lin_ref_table.loc[wildcards.lineage, 'File']),
        fasta = REF_FASTA,
        gff = REF_GFF
    output:
        REFDIR + "{lineage}_liftoff.gff_polished"
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
        "&> {log}"

rule ref2ref_agat:
    input: 
        fasta = REF_FASTA,
        ref_liftoff = REFDIR + "{lineage}_liftoff.gff_polished"
    output:
        cds = REFDIR + "{lineage}_predicted_cds.fa",
        prots = REFDIR + "{lineage}_predicted_proteins.fa"
    conda:
        "agat.yaml"
    log:
        cds = "logs/references/{lineage}_ref_agat_cds.log",   
        prots = "logs/references/{lineage}_ref_agat_prots.log"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.ref_liftoff} "
        "-f {input.fasta} "
        "-o {output.cds} "
        "&> {log.cds} "
        " && "
        "agat_sp_extract_sequences.pl "
        "-g {input.ref_liftoff} "
        "-f {input.fasta} "
        "-o {output.prots} "
        "-p &> {log.prots}"  
