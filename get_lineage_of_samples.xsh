import pandas as pd
from collections import defaultdict

# Get table with strain Name and SRS
sample_names = pd.read_table("PRJNA382844_samples.txt", header=None, names=["SAM", "Name", "SRS"])
sample_names.drop('SAM', axis = 1, inplace=True)
sample_names['Name'] = sample_names['Name'].str.replace("Cryptococcus neoformans var. grubii", "")
sample_names['Name'] = sample_names['Name'].str.replace("Cryptococcus_neoformans var. grubii_", "")
sample_names['Name'] = sample_names['Name'].str.replace("Cryptococcus_neoformans_", "")
sample_names['Name'] = sample_names['Name'].str.replace("Cryptococcus neoformans ", "")
sample_names['Name'] = sample_names['Name'].str.replace("Cryptococcus sp.", "")
sample_names.sort_values('Name', inplace=True)
sample_names['Name'] = sample_names['Name'].str.strip()

# Get table with name (Strain), lineage (Group) and SRA Accession
original_metadata = pd.read_csv("Desjardins_Supplemental_Table_S1.csv")
original_codes = original_metadata.drop(original_metadata.columns.difference(['Strain','Group', 'SRA Accession']), axis = 1)
original_codes.sort_values('Strain', inplace=True)

#Combine both tables and get samples with no coincidence in both
combined_table = sample_names.set_index('Name').join(original_codes.set_index('Strain'))
missing_samples = combined_table[combined_table['Group'].isna()]

# Get accession from missing SRSs
dfMissingSRS = pd.DataFrame([])
for srs in missing_samples.SRS:
        print(f"Processing {srs}")
        srr = $(esearch -db sra -query @(srs) | efetch -format runinfo -mode xml | xtract -pattern Row -element Sample Experiment Run SampleName)
        stringList = srr.split('\n')
        dfsrs = pd.DataFrame([item.split('\t') for item in stringList], columns = ['Sample', 'Experiment', 'Run', 'Name'])
        dfMissingSRS = pd.concat([dfMissingSRS, dfsrs])
dfMissingSRS.dropna(inplace=True)

#
original_codes[['SA1','SA2','SA3','SA4']] = original_codes['SRA Accession'].str.split(', ',expand=True)
original_codes.drop('SRA Accession', axis = 1, inplace = True)
original_codes = pd.lreshape(original_codes, {'SRA Accession': ['SA1', 'SA2', 'SA3', 'SA4']})

original_codes[['Run', 'Experiment']] = original_codes['SRA Accession'].str.split('SRX', expand=True)
original_codes[['Experiment']] = 'SRX'+original_codes[['Experiment']]
original_codes.drop('SRA Accession', axis = 1, inplace = True)

original_SRX = original_codes[original_codes['Experiment'].notna()]
original_SRX.drop('Run', axis=1,inplace=True)
combined_missing_SRX = dfMissingSRS.set_index('Experiment').join(original_SRX.set_index('Experiment'))
combined_missing_SRX = combined_missing_SRX[combined_missing_SRX['Group'].notna()]    

original_SRR = original_codes[original_codes['Run'] != '']
original_SRR.drop('Experiment', axis=1,inplace=True)
combined_missing_SRR = dfMissingSRS.set_index('Run').join(original_SRR.set_index('Run'))
combined_missing_SRR = combined_missing_SRR[combined_missing_SRR['Group'].notna()]    

complete_metadata = pd.concat([combined_missing_SRR, combined_missing_SRX])