configfile: "config.yaml"

import pandas as pd
samplefile=(pd.read_csv(config["sample_file"], sep=","))
samples=list(set(samplefile["sample"]))

#def get_sample_codes(wildcards):
#    return samples[wildcards.sample]

rule all:
    input:
        expand("fastq_combined/{sample}_1.fq.gz",sample=samples),
        expand("fastq_combined/{sample}_2.fq.gz",sample=samples)

rule combine_fastq:
    input:
        "test_read_pair_table.csv"
    output:
        "fastq_combined/{sample}_1.fq.gz",
        "fastq_combined/{sample}_2.fq.gz"
    shell:
        "xonsh fastq-combiner.xsh {{wildcard.sample}} {input} fastqs/ fastq_combined/" #Still does not work how to call the sample wildcard

