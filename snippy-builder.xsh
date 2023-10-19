# Snippy based pipeline

THREADS = 12
OUTDIR = "snippy-analysis"
SAMPLESTABLE = "sample_metadata.csv"
REFGENOMETABLE = "lineage_references.csv"
FASTQDIR = "fastq_combined"
REFDIR = "Reference_Genomes"

from collections import defaultdict
from pathlib import Path
import sys
import pandas as pd

# Create dataframe from SAMPLESTABLE and REFGENOMETABLE with the sample, read filenames, lineage and reference assembly file information

lineage = pd.read_csv(SAMPLESTABLE)
lineage = lineage[["Sample","Group"]]
lineage = lineage.rename(columns={'Sample': 'sample', 'Group': 'lineage'})

d={'sample': lineage["sample"], 'file1': lineage["sample"]+"_1.fq.gz", 'file2': lineage["sample"]+"_2.fq.gz"}
df = pd.DataFrame(data=d)
df = df.set_index('sample').join(lineage.set_index('sample'))
df.reset_index(inplace=True)

reference = pd.read_csv(REFGENOMETABLE)
reference = reference[["Group","File"]]
reference = reference.rename(columns={'Group': 'lineage', 'File': 'refgenome'})

df = df.set_index('lineage').join(reference.set_index('lineage'))
df.reset_index(inplace=True)


# Create a working directory where snippy results will be stored

workpath = Path(OUTDIR)
if not workpath.exists():
    workpath.mkdir()

# Run snippy for each sample taking the apropriate fq.gz files and reference genome according to lineage

for row in df.itertuples():
    sample, f1, f2 ,ref = row.sample, row.file1, row.file2, row.refgenome
    spath = workpath / sample
    if not workpath.exists():
        workpath.mkdir()
    if (spath / "snippy.done").exists():  # indicates already completed
        continue
    print(f"Analyzing {sample}")
    f1 = Path(FASTQDIR) / f1
    f2 = Path(FASTQDIR) / f2
    ref = Path(REFDIR) / ref
    if (not f1.exists()) or (not f2.exists()):
        print("One or more missing fastq files:  {f1}, {f2}")
        sys.exit(1)
    cmd = ["snippy",  
           f"--outdir={str(spath)}", 
           f"--cpus={THREADS}",  
           f"--ref={ref}", 
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

    


