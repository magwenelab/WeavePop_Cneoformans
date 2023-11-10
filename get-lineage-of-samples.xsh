#Inputs
METADATATABLE = "Desjardins_Supplemental_Table_S1.csv"
READSTABLE = "read_pair_table.csv"
#Outputs
OUTSAMPLESTABLE = "sample_metadata.csv"

import pandas as pd
import csv
from pathlib import Path  

# Get Run and Experiment accessions from SRS sample names from Entrez Direct
sample_names = pd.read_csv(READSTABLE)
sample_names = set(sample_names["sample"])
all_accessions = pd.DataFrame([])
for srs in sample_names:
        print(f"Processing {srs}")
        srr = $(esearch -db sra -query @(srs) | efetch -format runinfo -mode xml | xtract -pattern Row -element Sample Experiment Run SampleName)
        stringList = srr.split('\n')
        dfsrs = pd.DataFrame([item.split('\t') for item in stringList], columns = ['Sample', 'Experiment', 'Run', 'Name'])
        all_accessions = pd.concat([all_accessions, dfsrs])
all_accessions.dropna(inplace=True)

#Get origianl metadata from Desjardins table:
original_metadata = pd.read_csv(METADATATABLE)
original_codes = original_metadata.drop(original_metadata.columns.difference(['Strain','Group', 'SRA Accession']), axis = 1)
#Split SRA accession column because the same cell has many accessions and take to long format:
original_codes[['SA1','SA2','SA3','SA4']] = original_codes['SRA Accession'].str.split(', ',expand=True)
original_codes.drop('SRA Accession', axis = 1, inplace = True)
original_codes = pd.lreshape(original_codes, {'SRA Accession': ['SA1', 'SA2', 'SA3', 'SA4']})
#Separate accessions into Run and Experiment accessions:
original_codes[['Run', 'Experiment']] = original_codes['SRA Accession'].str.split('SRX', expand=True)
original_codes[['Experiment']] = 'SRX'+original_codes[['Experiment']]
original_codes.drop('SRA Accession', axis = 1, inplace = True)
#Combine Experiment accessions of the original table with the rest of the obtained accessions:
original_SRX = original_codes[original_codes['Experiment'].notna()]
original_SRX.drop('Run', axis=1,inplace=True)
combined_SRX = all_accessions.set_index('Experiment').join(original_SRX.set_index('Experiment'))
combined_SRX = combined_SRX[combined_SRX['Group'].notna()]
combined_SRX = combined_SRX.reset_index() 
#Combine Run accessions of the original table with the rest of the obtained accessions:
original_SRR = original_codes[original_codes['Run'] != '']
original_SRR.drop('Experiment', axis=1,inplace=True)
combined_SRR = all_accessions.set_index('Run').join(original_SRR.set_index('Run'))
combined_SRR = combined_SRR[combined_SRR['Group'].notna()]    
combined_SRR = combined_SRR.reset_index()
#Combine all accessions obtained from Runs and Experiments into one table:
accs_and_group = pd.concat([combined_SRR, combined_SRX])
#Maintain only SRS accession:
SRS_names = accs_and_group.drop(['Run','Experiment','Group'], axis =1)
SRS_names.drop_duplicates(inplace = True)
#Join all original metadata with SRS accessions:
sample_metadata = SRS_names.set_index('Strain').join(original_metadata.set_index('Strain'))
sample_metadata = sample_metadata.reset_index()
#Save all metadata with SRS names and names gotten from Entrez into csv file
filepath = Path(OUTSAMPLESTABLE)  
sample_metadata.to_csv(filepath, index=False)  