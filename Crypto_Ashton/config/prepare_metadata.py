#The paper_metadata.csv file was downloaded from https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-10092-5/MediaObjects/41467_2019_10092_MOESM5_ESM.xlsx
#It was converted to csv, the first two lines were removed and the values in Sequencing ID where used to fill the ENA accession column where it was missing
#The values in ENA accession of the Ashton Study were used to download the sequencing data with the downloading-tools workflow.
#The reads_tables.csv file was created by the downloading-tools workflow
#The run ERR2624135 was not possible to download.
import pandas as pd

original = pd.read_csv('config/paper_metadata.csv', header = 0)
original = original[original['ENA accession'] != 'ERR2624135']
reads = pd.read_csv('config/reads_table.csv', header = 0)
reads = reads[['sample', 'run']]
joined = pd.merge(original, reads, left_on = 'ENA accession', right_on = 'run', how = 'left')
joined = joined.drop(columns = ['run'])
joined = joined.rename(columns = {'Lab ID': 'strain'})
joined['group'] = 'VNI'
# Put the columns sample, group, strain to the left of the rest
joined = joined[['sample', 'group', 'strain'] + [col for col in joined.columns if col not in ['sample', 'group', 'strain']]]
# Filter out the Desjardins samples
joined = joined[joined['Study'] != 'Desjardins']
joined.to_csv('config/sample_metadata.csv', index = False, header = True)