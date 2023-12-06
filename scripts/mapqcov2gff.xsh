#/usr/bin/env xonsh -c
import pandas as pd
from pathlib import Path  
import click

@click.command()
@click.argument("mapqbed", type=click.Path(exists=True))
@click.argument("covbed", type=click.Path(exists=True))
@click.argument("gff", type=click.Path(exists=True))
@click.argument("covmapq", type=str)
@click.argument("newgff", type=str)

def mapqcov2gff(mapqbed, covbed, gff, covmapq, newgff):
    # Make merged version of BED file with MAPQ and Coverage columns
    mapq = pd.read_csv(mapqbed, names = ["Chromosome", "Start", "End", "MAPQ"], sep = "\t" )
    cov = pd.read_csv(covbed, names = ["Chromosome", "Start", "End", "Coverage"], sep = "\t" )
    df = pd.merge(mapq, cov, on=["Chromosome", "Start", "End"])
    filepath = Path(covmapq)  
    df.to_csv(filepath, index=False, header = False, sep='\t')  
    # Intersect MAPQ and Coverage of windows with GFF file
    gff = $(bedtools intersect -wa -wb -a @(gff)  -b @(covmapq))
    stringList = gff.split('\n')
    dfgff = pd.DataFrame([item.split('\t') for item in stringList], columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute", "Chromosome", "Start", "End", "MAPQ", "COV"])
    dfgff = dfgff.dropna()
    # Get average of MAPQ and Coverage of all windows covered by each feature
    dfgff.loc[:,'MAPQ'] = pd.to_numeric(dfgff['MAPQ'], errors='coerce')
    dfgff.loc[:,'COV'] = pd.to_numeric(dfgff['COV'], errors='coerce')
    dfgff['AvgMAPQ'] = dfgff['MAPQ'].groupby(dfgff['attribute']).transform('mean')
    dfgff['AvgCOV'] = dfgff['COV'].groupby(dfgff['attribute']).transform('mean')
    new_gff = dfgff.loc[:, ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute', 'AvgMAPQ', 'AvgCOV']]
    new_gff = new_gff.drop_duplicates()
    new_gff['AvgMAPQ'] = new_gff['AvgMAPQ'].apply(lambda x: round(x, 2)).astype(str)
    new_gff['AvgCOV'] = new_gff['AvgCOV'].apply(lambda x: round(x, 2)).astype(str)
    # Add average values to GFF and write to new GFF file
    new_gff['attribute'] = new_gff['attribute'] + ';window_avg_mapq=' + new_gff['AvgMAPQ']  + ';window_avg_cov=' + new_gff['AvgCOV']
    new_gff = new_gff.loc[:, ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']]
    filepath = Path(newgff)  
    new_gff.to_csv(filepath, index=False, header = False, sep='\t')

if __name__ == "__main__":
    mapqcov2gff() # Runs the function