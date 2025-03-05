#!/usr/bin/bash

# run this script over a set of directories with raw VCF files
# to get a quick summary of decent quality called het sites:
# tail -n +2 /FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Desjardins/config/metadata.csv | cut -d',' -f2 | while read line
# > do
# > bash ../../../scripts/het_stats.sh /FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Desjardins/results/01.Samples/snippy/$line/snps.raw.vcf $line
# > done
INFILE=$1
SAMPLE=$2
HETEROFILE=$SAMPLE.hetero.csv
HOMOFILE=$SAMPLE.homo.csv
SUMMARYFILE=summary.csv

echo chrom,pos,ref,alt,ref_freq > $HETEROFILE
cat $INFILE |
    bcftools view --min-alleles 2 --max-alleles 2 -v snps -g het -i 'FMT/DP>20 && QUAL>20 && INFO/AB > 0.2 && INFO/AB < 0.8' |
    bcftools query -f '%CHROM,%POS,%REF,%ALT,%INFO/AB\n' >> $HETEROFILE

echo chrom,pos,ref,alt,ref_freq > $HOMOFILE
cat $INFILE |
    bcftools view  -g ^het -i 'FMT/DP>20 && QUAL>20' |
    bcftools query -f '%CHROM,%POS,%REF,%ALT,%INFO/AB\n' >> $HOMOFILE

nhetero=$(tail -n +2 $HETEROFILE | wc -l )
nhomo=$(tail -n +2 $HOMOFILE | wc -l )
sum=$(($nhetero + $nhomo))
fhetero=$(($nhetero / $sum))
fhomo=$(($nhomo / $sum))

if ! [[ -s $SUMMARYFILE ]]; then
    echo "sample,n_hetero,n_homo,sum,freq_hetero,freq_homo" > $SUMMARYFILE
    echo $SAMPLE,$nhetero,$nhomo,$sum,$fhetero,$fhomo >> $SUMMARYFILE
else
    echo $SAMPLE,$nhetero,$nhomo,$sum,$fhetero,$fhomo >> $SUMMARYFILE

fi

