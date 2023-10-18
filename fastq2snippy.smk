configfile: "config.yaml"

import pandas as pd
sample_file = config["sample_file"]
SAMPLES = pd.read_table(sample_file, header = None)

rule all:
    input:
        expand("fastq_combined/{sample}_1.fq.gz",sample=SAMPLES),
        expand("fastq_combined/{sample}_2.fq.gz",sample=SAMPLES)

rule combine_fastq:
    input:
        expand("{sample}",sample=SAMPLES),
        csvfile = config["read_pair_table"] # Do something to put here the csv file needed for fastq-combiner.xsh
    output:
        "fastq_combined/{sample}_1.fq.gz"
        "fastq_combined/{sample}_2.fq.gz"
    run:
        "xonsh fastq-combiner.xsh {sample} {csvfile} fastqs/ fastq_combined/"