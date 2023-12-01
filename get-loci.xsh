#/usr/bin/env xonsh


# Get MAT loci protein_coding_gene IDs: Everything between SXI1 and STE12
# awk '/SXI1/,/STE12/' references/FungiDB-65_CneoformansH99.gff.tsv | grep protein_coding_gene | cut -f13 > results/MAT.txt
# Get rRNA IDs from the level2 primary_tag = rRNA and converting the ID into gene_ID (because the level1 primary tag is ncRNA_gene and the pattern rRNA is also in descriptions that are nor rRNA genes)
# awk '$3 == "rRNA" {print $0}' references/FungiDB-65_CneoformansH99.gff.tsv  | cut -f13 | cut -d'-' -f1 > results/rRNA.txt
# Get centromere delimiting-gene IDs from Janbon 2014 

import pandas as pd
from pathlib import Path

MAT = Path('results/MAT.txt')
CENTROMERE = Path('results/centromeres.txt')
RRNA = Path('results/rRNA.txt')
lineages = ["VNI", "VNII", "VNBI", "VNBII"]

mydfs = {}
for lin in lineages:
    mydfs[lin] = pd.read_csv(Path('references/' + lin + '_liftoff.gff_polished.tsv'), sep='\t', header=0)

annotations = pd.concat(mydfs)
genes = annotations[["seq_id", "primary_tag", "start", "description", "gene_id", "ID", "Name"]]
level1_bool = genes.primary_tag.str.contains("protein_coding_gene") | genes.primary_tag.str.contains("ncRNA_gene")
level1 = genes[level1_bool]

mat_list = list(pd.read_csv(MAT, header = None)[0])
centromere_list = list(pd.read_csv(CENTROMERE, header = None)[0])
rrna_list = list(pd.read_csv(rRNA, header = None)[0])

mat = level1[level1['ID'].isin(mat_list)]
rrna = level1[level1['ID'].isin(rrna_list)]
centromere = level1[level1['ID'].isin(centromere_list)]

mat['Loci'] = "MAT"
rrna['Loci'] = "rRNA"
centromere['Loci'] = "Centromere"

loci = [mat, rrna, centromere]
all_loci = pd.concat(loci)

all_loci.to_csv('results/loci_interest.tsv',index= False, sep ='\t')