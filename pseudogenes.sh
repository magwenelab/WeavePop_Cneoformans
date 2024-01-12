mkdir pseudogenes
head -n3 files/samples.txt | while read sample
do
    agat_sp_flag_premature_stop_codons.pl --gff analysis/${sample}/lifted.gff_polished --fasta analysis/${sample}/snps.consensus.fa -out pseudogenes/${sample}.gff 
    grep "contained" pseudogenes/${sample}_report.txt | cut -d " " -f5 > pseudogenes/${sample}_pseudogene_list.tsv
    sed -i "s/$/,\\${sample}/" pseudogenes/${sample}_pseudogene_list.tsv
done

cat pseudogenes/*.tsv > pseudogenes_list.tsv


