#/usr/bin/env xonsh

for d in p"Samples/".resolve().glob("SRS*/SRR*"):
    if (not len(list(d.glob("*.sra")))) and (not len(list(d.glob("*.vdbcache")))):
        print(f"Retrieving: {d.parent.name}/{d.name}")
        prefetch @(d.name) -O @(str(d.parent))
