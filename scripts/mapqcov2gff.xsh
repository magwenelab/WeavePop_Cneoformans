#/usr/bin/env xonsh -c
import pandas as pd
from pathlib import Path  

MAPQBED = "mapq_window.bed"
COVBED = "coverage.regions.bed"
COVMAPQ = "mapq_cov_window.bed"
GFF = "lifted.gff_polished"
NEWGFF = "annotation.gff"

mapq = pd.read_csv(MAPQBED, names = ["Chromosome", "Start", "End", "MAPQ"], sep = "\t" )
cov = pd.read_csv(COVBED, names = ["Chromosome", "Start", "End", "Coverage"], sep = "\t" )

df = pd.merge(mapq, cov, on=["Chromosome", "Start", "End"])

filepath = Path(COVMAPQ)  
df.to_csv(filepath, index=False, header = False, sep='\t')  

gff = $(bedtools intersect -wa -wb -a @(GFF)  -b @(COVMAPQ))
stringList = gff.split('\n')
dfgff = pd.DataFrame([item.split('\t') for item in stringList], columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute", "Chromosome", "Start", "End", "MAPQ", "COV"])

dfgff = dfgff.dropna()

dfgff.loc[:,'MAPQ'] = pd.to_numeric(dfgff['MAPQ'], errors='coerce').round(2)
dfgff.loc[:,'COV'] = pd.to_numeric(dfgff['COV'], errors='coerce')
dfgff['AvgMAPQ'] = dfgff['MAPQ'].groupby(dfgff['attribute']).transform('mean')
dfgff['AvgCOV'] = dfgff['COV'].groupby(dfgff['attribute']).transform('mean')

new_gff = dfgff.loc[:, ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute', 'AvgMAPQ', 'AvgCOV']]
new_gff = new_gff.drop_duplicates()

new_gff['AvgMAPQ'] = new_gff['AvgMAPQ'].apply(lambda x: round(x, 2)).astype(str)
new_gff['AvgCOV'] = new_gff['AvgCOV'].apply(lambda x: round(x, 2)).astype(str)


new_gff['attribute'] = new_gff['attribute'] + ';window_avg_mapq=' + new_gff['AvgMAPQ']  + ';window_avg_cov=' + new_gff['AvgCOV']
new_gff = new_gff.loc[:, ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']]

filepath = Path(NEWGFF)  
new_gff.to_csv(filepath, index=False, header = False, sep='\t')