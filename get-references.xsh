#Inputs
SAMPLESTABLE = "sample_metadata.csv"
REFGENOMETABLE = "lineage_references.csv"
#Outputs
OUTSAMPLEREFERENCE = "sample_reference.csv"

from collections import defaultdict
from pathlib import Path
import sys
import pandas as pd

# Create dataframe from SAMPLESTABLE and REFGENOMETABLE with the sample, read filenames, lineage and reference assembly file information

lineage = pd.read_csv(SAMPLESTABLE)
lineage = lineage[["sample","group"]]
lineage = lineage.rename(columns={'group': 'lineage'})

d={'sample': lineage["sample"], 'file1': lineage["sample"]+"_1.fq.gz", 'file2': lineage["sample"]+"_2.fq.gz"}
df = pd.DataFrame(data=d)
df = df.set_index('sample').join(lineage.set_index('sample'))
df.reset_index(inplace=True)

reference = pd.read_csv(REFGENOMETABLE)
reference = reference[["group","File"]]
reference = reference.rename(columns={'group': 'lineage', 'File': 'refgenome'})

df = df.set_index('lineage').join(reference.set_index('lineage'))
df.reset_index(inplace=True)

filepath = Path(OUTsampleREFERENCE)  
df.to_csv(filepath, index=False)