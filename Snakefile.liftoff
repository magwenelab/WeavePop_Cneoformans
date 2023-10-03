from pathlib import Path
SAMPLES = [str(p) for p in Path().glob("SRS*")]

REF_FASTA = "FungiDB-53_CneoformansH99_Genome.fasta"
REF_GFF = "FungiDB-53_CneoformansH99_PMM.gff"

rule all:
    input:
        expand("{sample}/lifted.gff_polished",sample=SAMPLES),
        # expand("{sample}/predicted_cds.fa", sample=SAMPLES),
        # expand("{sample}/predicted_proteins.fa", sample=SAMPLES),
        # expand("{sample}/predicted_proteins.fa.fai", sample=SAMPLES),
        # "reference_predicted_proteins.fa"  

rule by_cds:
    input:
        prots = "protein_list.txt",
        samples = "samples.txt"
    output:
        "cds/done.txt"
    run:
        with open(input.prots, 'r') as f:
            for line in f:
                prot = line.strip()
                shell(f"./getcds.sh {prot}")
        shell(f'touch {output}')    

rule by_protein:
    input:
        prots = "protein_list.txt",
        samples = "samples.txt"
    output:
        "proteins/done.txt"
    run:
        with open(input.prots, 'r') as f:
            for line in f:
                prot = line.strip()
                shell(f"./getprots.sh {prot}")
        shell(f'touch {output}')                


rule index_proteins:
    input:
        "{sample}/predicted_proteins.fa"
    output:
        "{sample}/predicted_proteins.fa.fai"
    shell:
        "seqkit faidx {input}"


rule sample_list:
    output:
        "samples.txt"
    shell:
        "ls -d SRS* > {output}"

rule protein_list:
    input:
        "reference_predicted_proteins.fa"
    output:
        "protein_list.txt"
    shell:
        "seqkit seq -n -i {input} > {output}"

rule ref2ref_agat:
    input: 
        fasta = REF_FASTA,
        ref_liftoff = "ref_liftoff.gff_polished"
    output:
        cds = "reference_predicted_cds.fa",
        aa = "reference_predicted_proteins.fa"        
    conda:
        "agat.yaml"   
    shell:
        """agat_sp_extract_sequences.pl \
            -g {input.ref_liftoff} \
            -f {input.fasta}         \
            -o {output.cds}           \
        && \       
        agat_sp_extract_sequences.pl \
            -g ref_liftoff.gff_polished  \
            -f {input.fasta}         \
            -o {output.aa} -p""" 


rule ref2ref_liftoff:
    input:
        fasta = REF_FASTA,
        gff = REF_GFF
    output:
        "ref_liftoff.gff_polished"
    shell:
        """liftoff            \
            -g {input.gff}  \
            -polish         \
            -o ref_liftoff.gff  \
            {input.fasta}   \
            {input.fasta}""" 
            


rule agat:
    input:
        "{sample}/liftoff.done",
        "{sample}/lifted.gff_polished"
    output:
        "{sample}/predicted_cds.fa",
        "{sample}/predicted_proteins.fa"
    conda:
        "agat.yaml"
    shell:
        "seqkit seq -w 70 {wildcards.sample}/snps.consensus.fa > "
        "{wildcards.sample}/snps.consensus.wrapped.fa && "
        "agat_sp_extract_sequences.pl "
        "-g {input[1]} " 
        "-f {wildcards.sample}/snps.consensus.wrapped.fa "
        "-o {output[0]} "
        " && "
        "agat_sp_extract_sequences.pl "
        "-g {input[1]} " 
        "-f {wildcards.sample}/snps.consensus.wrapped.fa "
        "-o {output[1]} "
        "-p "        


rule liftoff:
    input:
        "{sample}/snps.consensus.fa"
    params:
        refgenome = "{sample}/ref.fa",
        refgff = f"{REF_GFF}"
    output:
        "{sample}/lifted.gff",        
        "{sample}/lifted.gff_polished",
        "{sample}/liftoff.done"
    shell:
        "cp {params.refgff} {wildcards.sample}/ref.gff && "
        "liftoff "
        "-g {wildcards.sample}/ref.gff "
        "-polish "
        "-dir {wildcards.sample}/liftoff "
        "-o {output[0]} "
        "{input} "
        "{params.refgenome} "
        " && touch {output[2]} "


