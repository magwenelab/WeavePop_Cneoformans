configfile: "config.yaml"

import pandas as pd
samplefile=(pd.read_csv(config["sample_file"], sep=","))
samples=list(set(samplefile["sample"]))

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
