import pandas as pd
from collections import defaultdict

sample_names = pd.read_table("PRJNA382844_samples.txt", header=None, names=["SAM", "Name", "SRS"])
sample_names.drop('SAM', axis = 1, inplace=True)

sample_names['Name'] = sample_names['Name'].str.replace("Cryptococcus neoformans var. grubii", "")
sample_names['Name'] = sample_names['Name'].str.replace("Cryptococcus_neoformans var. grubii_", "")
sample_names['Name'] = sample_names['Name'].str.replace("Cryptococcus_neoformans_", "")
sample_names['Name'] = sample_names['Name'].str.replace("Cryptococcus neoformans ", "")
sample_names['Name'] = sample_names['Name'].str.replace("Cryptococcus sp.", "")

original_metadata = pd.read_csv("Desjardins_Supplemental_Table_S1.csv")
original_metadata.drop(original_metadata.columns.difference(['Strain','Group', 'SRA Accession']), axis = 1, inplace=True)

sample_names.sort_values('Name', inplace=True)
original_metadata.sort_values('Strain', inplace=True)

sample_names['Name'].equals(original_metadata['Strain'])

sample_names['Name'] = sample_names['Name'].str.strip()

combined_table = sample_names.set_index('Name').join(original_metadata.set_index('Strain'))

missing_samples = combined_table[combined_table['Group'].isna()]

dfMissingSRS = pd.DataFrame([])
for srs in missing_samples.SRS:
        print(f"Processing {srs}")
        srr = $(esearch -db sra -query @(srs) | efetch -format runinfo -mode xml | xtract -pattern Row -element Sample Experiment Run SampleName)
        stringList = srr.split('\n')
        dfsrs = pd.DataFrame([item.split('\t') for item in stringList], columns = ['Sample', 'Experiment', 'Run', 'Name'])
        dfMissingSRS = pd.concat([dfMissingSRS, dfsrs])
dfMissingSRS.dropna(inplace=True)

original_metadata[['SA1','SA2','SA3','SA4']] = original_metadata['SRA Accession'].str.split(', ',expand=True)
original_metadata.drop('SRA Accession', axis = 1, inplace = True)
df1 = pd.lreshape(original_metadata, {'SRA Accession': ['SA1', 'SA2', 'SA3', 'SA4']})

