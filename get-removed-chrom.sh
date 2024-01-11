#/bin/bash
#
# Remove chromosome/contig from fasta file
# Keeping the original as name.fasta.original

scriptargs="path-to-fasta seq_id"
nexpargs=2
nargs=$#

if [[ $nargs -ne $nexpargs ]]
then
    echo
    echo "Usage: $(basename $0) $scriptargs"
    echo
    exit
fi

fasta=$1
seqid=$2

echo "Removing $seqid from $fasta"
seqkit grep -n -r -v -p ${seqid} ${fasta} > ${fasta}.modif
mv ${fasta} ${fasta}.original
mv ${fasta}.modif ${fasta}
