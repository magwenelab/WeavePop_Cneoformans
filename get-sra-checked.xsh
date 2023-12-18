#/usr/bin/env xonsh

for d in p"srafiles/".resolve().glob("SRS*/SRR*"):
    if (not len(list(d.glob("*.sra")))) and (not len(list(d.glob("*.vdbcache")))):
        print(f"Prefetching: {d.parent.name}/{d.name}")
        prefetch @(d.name) -O @(str(d.parent))
