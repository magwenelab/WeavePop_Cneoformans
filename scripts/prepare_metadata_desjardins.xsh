# Run in environment `sra-tools` with the following commands:
# conda activate sra-tools
# xonsh scripts/prepare_metadata_desjardins.xsh

import pandas as pd
import csv
from pathlib import Path  
import numpy as np

#Inputs
METADATATABLE = "Crypto_Desjardins/config/Desjardins_Supplemental_Table_S1.csv"
READSTABLE = "Crypto_Desjardins/config/reads_table.csv"
#Outputs
OUTSAMPLESTABLE = "Crypto_Desjardins/config/metadata.csv"

# Get Run and Experiment accessions from SRS sample names from Entrez Direct
sample_names = pd.read_csv(READSTABLE)
sample_names = set(sample_names["sample"])
all_accessions = pd.DataFrame([])
for srs in sample_names:
        print(f"Processing {srs}")
        attempt = 0
        srr = ''
        while attempt < 20:
                srr = $(esearch -db sra -query @(srs) | efetch -format runinfo -mode xml | xtract -pattern Row -element Sample Experiment Run SampleName)
                if srr.strip():
                        break
                print(f"Attempt {attempt+1} failed for {srs}, retrying...")
                attempt += 1
        if not srr.strip():
                print(f"Failed to retrieve data for {srs} after 20 attempts. Skipping.")
                continue
        stringList = srr.split('\n')
        dfsrs = pd.DataFrame([item.split('\t') for item in stringList], columns = ['Sample', 'Experiment', 'Run', 'Name'])
        all_accessions = pd.concat([all_accessions, dfsrs])
all_accessions.dropna(inplace=True)

#Get origianl metadata from Desjardins table:
original_metadata = pd.read_csv(METADATATABLE)
original_codes = original_metadata.drop(original_metadata.columns.difference(['Strain','Group', 'SRA Accession']), axis = 1)

#Split SRA accession column because the same cell has many accessions and convert to long format:
original_codes[['SA1','SA2','SA3','SA4']] = original_codes['SRA Accession'].str.split(', ',expand=True)
original_codes.drop('SRA Accession', axis = 1, inplace = True)
original_codes = pd.lreshape(original_codes, {'SRA Accession': ['SA1', 'SA2', 'SA3', 'SA4']})

#Separate accessions into Run and Experiment accessions:
original_codes[['Run', 'Experiment']] = original_codes['SRA Accession'].str.split('SRX', expand=True)
original_codes[['Experiment']] = 'SRX'+original_codes[['Experiment']]
original_codes.drop('SRA Accession', axis = 1, inplace = True)

#Combine Experiment accessions of the original table with the rest of the obtained accessions:
original_SRX = original_codes[original_codes['Experiment'].notna()].copy()
original_SRX.drop('Run', axis=1,inplace=True)
combined_SRX = all_accessions.set_index('Experiment').join(original_SRX.set_index('Experiment'))
combined_SRX = combined_SRX[combined_SRX['Group'].notna()]
combined_SRX = combined_SRX.reset_index() 

#Combine Run accessions of the original table with the rest of the obtained accessions:
original_SRR = original_codes[original_codes['Run'] != ''].copy()
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

sample_metadata.rename(columns={"Sample": "sample",
                                "Strain": "strain",
                                "Group": "lineage",
                                "VNI subdivision": "vni_subdivision",
                                "Mating type": "mating_type", 
                                "Country of origin": "country_of_origin",
                                "Name": "name",
                                "Isolation source": "isolation_source",
                                "Broad Project": "broad_project",
                                "SRA Accession": "sra_accession", 
                                "Strain description": "strain_description",
                                "In GWAS":"in_gwas",
                                "Phenotyped": "phenotyped"}, inplace=True)

# Create source column based on isolation_source
sample_metadata['source'] = np.where(sample_metadata.isolation_source.str.contains("tree",case=False), "Environmental",
                   np.where(sample_metadata.isolation_source.str.contains("avian",case=False), "Environmental",
                   np.where(sample_metadata.isolation_source.str.contains("soil",case=False), "Environmental",
                   np.where(sample_metadata.isolation_source.str.contains("guano",case=False), "Environmental", "Clinical"))))

# Remove leading and trailing whitespace from all columns
sample_metadata = sample_metadata.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

# Replace not assigned values with NaN in mating_type column
sample_metadata['mating_type'] = sample_metadata['mating_type'].replace({'not assigned': np.nan})

# Reorder columns
sample_metadata = sample_metadata[['sample', 'strain', 'lineage', 'vni_subdivision', 'source', 'mating_type', 
                                   'country_of_origin', 'name', 'isolation_source', 'broad_project', 
                                   'sra_accession', 'strain_description', 'in_gwas', 'phenotyped']]

#Save all metadata with SRS names and names gotten from Entrez into csv file
filepath = Path(OUTSAMPLESTABLE)  
sample_metadata.to_csv(filepath, index=False)  