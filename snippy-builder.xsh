# Snippy based pipeline

THREADS = 12
OUTDIR = "snippy-analysis"
READPAIRTABLE = "largest_read_pair_table.csv"
REFGENOME = "Reference_Genomes/FungiDB-53_CneoformansH99_Genome.fasta"
FASTQDIR = "fastqs"

from collections import defaultdict
from pathlib import Path
import sys

import pandas as pd

# Read sample -> read pairs table from CSV file
df = pd.read_csv(READPAIRTABLE)

# Create a working directory where snippy results will be stored
# - dir name: snippy-analysis

workpath = Path(OUTDIR)
if not workpath.exists():
    workpath.mkdir()

for row in df.itertuples():
    sample, f1, f2 = row.sample, row.file1, row.file2
    spath = workpath / sample
    if not workpath.exists():
        workpath.mkdir()
    if (spath / "snippy.done").exists():  # indicates already completed
        continue
    print(f"Analyzing {sample}")
    f1 = Path(FASTQDIR) / f1
    f2 = Path(FASTQDIR) / f2
    if (not f1.exists()) or (not f2.exists()):
        print("One or more missing fastq files:  {f1}, {f2}")
        sys.exit(1)
    cmd = ["snippy",  
           f"--outdir={str(spath)}", 
           f"--cpus={THREADS}",  
           f"--ref={REFGENOME}", 
           f"--R1={f1}",
           f"--R2={f2}",
    ]

    result = !(@(cmd))
    if result:
        touch @(str(spath / "snippy.done"))  # write file to indicate succesful completion



# # For every Sample directory (SRS* directories) under "Samples/"
# # - get Sample name and path
# # - Find all records (SRR* files) in that Sample path
# #     - Sort the records by alphanumerical order
# #     - Pick the first records as the representative sequence record to work with

# Samples = p"Samples"
# sample_paths = list(Samples.glob("SRS*"))
# sample_names = [p.name for p in sample_paths]

# records = []
# sample2record = {}

# for p in sample_paths:
#     allrecs = list(p.resolve().glob("SRR*"))
#     rec = sorted(allrecs)[0]
#     records.append(rec)
#     sample2record[p.name] = rec.name
 
# # For each representative sequence record
# # - create a directory /snippy-analysis/SRS*/
# # - get the SRA file (*.sra) path
# # - create fastq files from that path using fasterq-dump (part of SRA tools), writing files to `snippy-analysis/SRS*/SRR*.sra_{1,2}.fastq`

# for (sample, run) in sample2record.items():
#     fastqdir = workpath / sample
#     if not fastqdir.exists():
#         fastqdir.mkdir()

#     srapath = Samples / sample / run / f"{run}.sra"
    
#     if not srapath.exists():
#         print(f"WARNING: Missing .sra file: {srapath}")
#         continue
    
#     print(f"Processing {srapath}:")
#     fasterq-dump @(srapath) -O @(fastqdir) -p

    


