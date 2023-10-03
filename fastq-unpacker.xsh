#/usr/bin/env xonsh

for each in p"Samples".resolve().rglob("*.sra"):
    ck = p"fastqs" / f"{each.parent.name}_1.fastq"
    if ck.exists(): 
        continue
    print(f"Processing {each}")
    fasterq-dump @(str(each)) -O "fastqs" -p -S