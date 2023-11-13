#/usr/bin/env xonsh

import pandas as pd
import csv
from pathlib import Path  

LINEAGETABLE = "./lineage_references.csv"

lineages = pd.read_csv(LINEAGETABLE)

paths = lineages['File'].tolist() 
