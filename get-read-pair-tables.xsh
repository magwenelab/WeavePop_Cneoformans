#/usr/bin/env xonsh

import csv

paired = []
unpaired = []

for srs in sorted(p"Samples".glob("SRS*")):
    for srr in sorted(srs.glob("SRR*")):
        r1 = p"fastqs" / f"{srr.name}_1.fastq"
        r2 = p"fastqs" / f"{srr.name}_2.fastq"
        if r1.exists() and r2.exists(): 
            size1 = r1.stat().st_size
            size2 = r2.stat().st_size
            paired.append((srs.name, srr.name, r1.name, r2.name, size1+size2))
        else:
            ufiles = ",".join([str(i.name) for i in p"fastqs".glob(f"{srr.name}*.fastq")])
            unpaired.append((srs.name, srr.name, ufiles))

with open('read_pair_table.csv', 'w', newline='') as csvfile:
    w = csv.writer(csvfile)
    w.writerow(("sample", "run", "file1", "file2", "size"))
    w.writerows(paired)


with open('unpaired_fastqs.csv', 'w', newline='') as csvfile:
    w = csv.writer(csvfile)
    w.writerow(("sample", "run", "files"))
    w.writerows(unpaired)