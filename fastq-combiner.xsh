#/usr/bin/env xonsh -c

from pathlib import Path
import pandas as pd
import click


@click.command()
@click.argument("sample", type=str)
@click.argument("csvfile", type=click.Path(exists=True))
@click.argument("inputdir", type=click.Path(exists=True))
@click.argument("outputdir")

def combine(sample, csvfile, inputdir, outputdir):
    """This script combines all fastq _1 and _2 files of one sample into only one _1.fq.gz and one _2.fq.gz.
    
    SAMPLE is the sample name.
    
    CSVFILE is a csv table with the fields "sample", "file1" and "file2" with the sample name, _1.fastq and _2.fastq filenames, respectively.
    
    INPUTDIR is the name of the directory where the _1.fastq and _2.fastq files are saved. Not compressed!
    
    OUTPUTDIR is the name of the new directory for the combined .fq.gz files. It will be created if necessary.
    """
    workpath = Path(outputdir)
    if not workpath.exists():
        workpath.mkdir()
    df = pd.read_csv(csvfile)
    sampledf = df[df["sample"] == sample]
    pinput, poutput = Path(inputdir), Path(outputdir)
    files1 = [str(pinput / i) for i in sampledf["file1"]]
    files2 = [str(pinput / i) for i in sampledf["file2"]]
    echo "Combining sample" @(sample)
    cat @(files1) | gzip -1 > @(str(poutput / f"{sample}_1.fq.gz"))
    cat @(files2) | gzip -1 > @(str(poutput / f"{sample}_2.fq.gz"))


if __name__ == "__main__":
    combine()