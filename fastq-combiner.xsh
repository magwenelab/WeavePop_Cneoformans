#/usr/bin/env xonsh -c

from pathlib import Path
import pandas as pd
import click


@click.command()
@click.argument("sample", type=str)
@click.argument("csvfile", type=click.Path(exists=True))
@click.argument("inputdir", type=click.Path(exists=True))
@click.argument("outputdir", type=click.Path(exists=True))
def combine(sample, csvfile, inputdir, outputdir):
    df = pd.read_csv(csvfile)
    sampledf = df[df["sample"] == sample]
    pinput, poutput = Path(inputdir), Path(outputdir)
    files1 = [str(pinput / i) for i in sampledf["file1"]]
    files2 = [str(pinput / i) for i in sampledf["file2"]]
    cat @(files1) | gzip -1 > @(str(poutput / f"{sample}_1.fq.gz"))
    cat @(files2) | gzip -1 > @(str(poutput / f"{sample}_2.fq.gz"))


if __name__ == "__main__":
    combine()