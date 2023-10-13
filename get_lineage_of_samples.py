import pandas as pd

sample_names = pd.read_table("PRJNA382844_samples.txt", header=None, names=["SAM", "Name", "SRS"])
sample_names.drop('SAM', axis = 1, inplace=True)

sample_names['Name'] = sample_names['Name'].str.replace("Cryptococcus neoformans var. grubii", "")
sample_names['Name'] = sample_names['Name'].str.replace("Cryptococcus_neoformans var. grubii_", "")
sample_names['Name'] = sample_names['Name'].str.replace("Cryptococcus_neoformans_", "")
sample_names['Name'] = sample_names['Name'].str.replace("Cryptococcus neoformans ", "")
sample_names['Name'] = sample_names['Name'].str.replace("Cryptococcus sp.", "")

original_metadata = pd.read_csv("Desjardins_Supplemental_Table_S1.csv")
original_metadata.drop(original_metadata.columns.difference(['Strain','Group']), axis = 1, inplace=True)

sample_names.sort_values('Name', inplace=True)
original_metadata.sort_values('Strain', inplace=True)

sample_names['Name'].equals(original_metadata['Strain'])

sample_names['Name'] = sample_names['Name'].str.strip()


