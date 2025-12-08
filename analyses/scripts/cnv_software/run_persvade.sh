#!/bin/bash

cd /FastData/czirion/WeavePop_Cneoformans/analyses/results/cnv_software
mkdir -p results_HMMcopy_AneuFinder/logs_HMMcopy_AneuFinder

singularity exec \
    -B /FastData/czirion/WeavePop_Cneoformans/ \
    -e ./mikischikora_persvade_v1.02.6.sif \
    bash -c \
    "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate perSVade_env && \
    parallel -j 30 --colsep '\t' -a /FastData/czirion/WeavePop_Cneoformans/analyses/results/tables/haploid_samples_with_cnv.tsv\
    'bash /FastData/czirion/WeavePop_Cneoformans/analyses/results/cnv_software/scripts/persvade_commands.sh {1} {3} &> results_HMMcopy_AneuFinder/logs_HMMcopy_AneuFinder/{1}.log'"
