configfile: "/analysis/czirion/CryptoDiversity/config.yaml"

import pandas as pd
samplefile=(pd.read_csv(config["sample_file"], sep=","))
SAMPS=list(set(samplefile["sample"]))

ref_table = (pd.read_csv(config["sample_reference_file"], sep=","))
ref_table.set_index('sample', inplace=True)

lin_ref_table = (pd.read_csv(config["lineage_reference_file"], sep=","))
LINS=list(lin_ref_table["Group"])
lin_ref_table.set_index('Group', inplace=True)

REF_FASTA = "FungiDB-53_CneoformansH99_Genome.fasta"
REF_GFF = "FungiDB-53_CneoformansH99_PMM.gff"

rule all:
    input:
        expand("{lineage}_liftoff.gff_polished",lineage=LINS),
        expand("{lineage}_predicted_proteins.fa",lineage=LINS),
        expand("{lineage}_predicted_cds.fa",lineage=LINS),
        expand("{sample}/liftoff/lifted.gff_polished", sample=SAMPS)

rule ref2ref_liftoff:
    input:
        target_refs = lambda wildcards: lin_ref_table.loc[wildcards.lineage, 'File'],
        fasta = REF_FASTA,
        gff = REF_GFF
    output:
        "{lineage}_liftoff.gff_polished"
    log:
        "logs/liftoff/ref_{lineage}.log" 
    shell:
        "liftoff "
        "-g {input.gff} "
        "-polish "
        "-o {wildcards.lineage}_liftoff.gff "
        "{input.target_refs} {input.fasta} "
        "&> {log}"

rule ref2ref_agat_cds:
    input: 
        fasta = REF_FASTA,
        ref_liftoff = "{lineage}_liftoff.gff_polished"
    output:
        cds = "{lineage}_predicted_cds.fa"      
    conda:
        "/home/czirion/miniconda3/envs/agat"   
    shell:
        """agat_sp_extract_sequences.pl \
            -g {input.ref_liftoff} \
            -f {input.fasta}         \
            -o {output.cds}"""

rule ref2ref_agat_prot:
    input: 
        fasta = REF_FASTA,
        ref_liftoff = "{lineage}_liftoff.gff_polished"
    output:
        aa = "{lineage}_predicted_proteins.fa"        
    conda:
        "/home/czirion/miniconda3/envs/agat"   
    shell:
        """agat_sp_extract_sequences.pl \
            -g {wildcards.lineage}_liftoff.gff_polished  \
            -f {input.fasta}         \
            -o {output.aa} -p"""   


rule liftoff:
    input:
        "{sample}/snps.consensus.fa"
    params:
        lin_ref = lambda wildcards: ref_table.loc[wildcards.sample, 'lineage'],
        refgenome = "{sample}/ref.fa"
        #refgff = "{lin_ref}_liftoff.gff_polished"
    output:
        "{sample}/liftoff/lifted.gff",        
        "{sample}/liftoff/lifted.gff_polished"
        #"{sample}/liftoff/liftoff.done"
    shell:
        #"cp {params.refgff} {wildcards.sample}/ref.gff && "
        "liftoff "
        "-g {params.lin_ref}_liftoff.gff_polished "
        "-polish "
        "-dir {wildcards.sample}/liftoff "
        "-o {wildcards.sample}/liftoff/lifted.gff "
        "{input} "
        "{params.refgenome} "                    
