configfile: "/analysis/czirion/CryptoDiversity/config.yaml"

import pandas as pd

lin_ref_table = (pd.read_csv(config["lineage_reference_file"], sep=","))
LINS=list(lin_ref_table["Group"])
lin_ref_table.set_index('Group', inplace=True)

REF_FASTA = "Reference_Genomes/FungiDB-53_CneoformansH99_Genome.fasta"
REF_GFF = "Reference_Genomes/FungiDB-53_CneoformansH99_PMM.gff"

rule all:
    input:
        expand("Reference_Genomes/{lineage}_liftoff.gff_polished",lineage=LINS),
        expand("Reference_Genomes/{lineage}_predicted_proteins.fa",lineage=LINS),
        expand("Reference_Genomes/{lineage}_predicted_cds.fa",lineage=LINS)

rule ref2ref_liftoff:
    input:
        target_refs = lambda wildcards: ("Reference_Genomes/" +  lin_ref_table.loc[wildcards.lineage, 'File']),
        fasta = REF_FASTA,
        gff = REF_GFF
    output:
        "Reference_Genomes/{lineage}_liftoff.gff_polished"
    log:
        "logs/{lineage}_ref_liftoff.log" 
    shell:
        "liftoff "
        "-g {input.gff} "
        "-polish "
        "-o Reference_Genomes/{wildcards.lineage}_liftoff.gff "
        "-dir Reference_Genomes/{wildcards.lineage}_intermediate_files "
        "-u Reference_Genomes/{wildcards.lineage}_unmapped_features.txt "
        "{input.target_refs} {input.fasta} "
        "&> {log}"

rule ref2ref_agat_cds:
    input: 
        fasta = REF_FASTA,
        ref_liftoff = "Reference_Genomes/{lineage}_liftoff.gff_polished"
    output:
        cds = "Reference_Genomes/{lineage}_predicted_cds.fa"    
    conda:
        "agat.yaml"  
    log:
        "logs/{lineage}_ref_agat_cds.log"   
    shell:
        """agat_sp_extract_sequences.pl \
            -g {input.ref_liftoff} \
            -f {input.fasta}         \
            -o {output.cds} &> {log}"""

rule ref2ref_agat_prot:
    input: 
        fasta = REF_FASTA,
        ref_liftoff = "Reference_Genomes/{lineage}_liftoff.gff_polished"
    output:
        aa = "Reference_Genomes/{lineage}_predicted_proteins.fa"        
    conda:
        "agat.yaml" 
    log:
        "logs/{lineage}_ref_agat_prots.log"      
    shell:
        """agat_sp_extract_sequences.pl \
            -g Reference_Genomes/{wildcards.lineage}_liftoff.gff_polished  \
            -f {input.fasta}         \
            -o {output.aa} -p &> {log}"""   
