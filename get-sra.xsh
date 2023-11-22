#/usr/bin/env xonsh
"""
Based on code examples here:
    https://github.com/schultzm/entrez_direct_tut
"""

from pathlib import Path
from collections import defaultdict

import pandas as pd

# BioProject identifier
project_id = "PRJNA382844"

sample_file = f"{project_id}_samples.txt"

# Get Accession numbers, Names, and SRS numbers for each sample
# associated with the bioproject and write this info to a file
print("Getting Sample Info")
if not Path(sample_file).exists():
    esearch -db bioproject -query @(project_id) | elink -target biosample | efetch -format docsum | xtract -pattern DocumentSummary -block Ids -element Id -group sra > @(sample_file)

# Read the Sample file as a data frame
df = pd.read_table(sample_file, header=None, names=["SAM", "Name", "SRS"])


# Get sequencing run IDs associated with each sample
# Don't reconstruct SRS to SRR table if it already exists
srs_file = f"{project_id}_SRStoSRR.txt"
dsrs = defaultdict(list)

print("Getting Run Info")
if not Path(srs_file).exists():
    ofile = open(srs_file,"w")
    for srs in df.SRS:
        print(f"Processing {srs}")
        srr = $(esearch -db sra -query @(srs) | efetch -format runinfo -mode xml | xtract -pattern Run -element Run)
        for item in srr.split("\n"):
            if not len(item): 
                continue
            dsrs[srs].append(item)
            ofile.write(f"{srs}\t{item}\n")
        ofile.flush()
    ofile.close()
else:
    dfsrs = pd.read_table(srs_file, header=None, names = ["SRS","SRR"])
    for (srs,srr) in zip(dfsrs.SRS, dfsrs.SRR):
        dsrs[srs].append(srr)


# Prefetch the SRA files
for sample, run in dsrs.items():
    mkdir -p Samples/@(sample)
    cd Samples/@(sample)
    prefetch @(run)
    cd ../..

