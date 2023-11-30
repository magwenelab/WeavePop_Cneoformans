#VNI
cut -f1,3,4,5,7,11,15,16,21 references/VNI_liftoff.gff_polished.tsv | head -n1 > results/VNI_genes.tsv
cut -f1,3,4,5,7,11,15,16,21 references/VNI_liftoff.gff_polished.tsv | grep -e protein_coding_gene -e ncRNA_gene -e pseudogene >> results/VNI_genes.tsv
 
head -n1 results/VNI_genes.tsv > results/loci_interest.tsv
#sed -i "s/$/\tLoci/" results/loci_interest.tsv

awk '/SXI1/,/STE12/' results/VNI_genes.tsv >> results/loci_interest.tsv
#sed -i "s/$/\tMAT/" results/loci_interest.tsv

cut -d',' -f2 results/centromeres.csv | while read line
do
awk -v i="$line" '$0 ~ i {print $0}' results/VNI_genes.tsv
done >> results/loci_interest.tsv

#VNBI
cut -f1,3,4,5,7,11,15,16,22 references/VNBI_liftoff.gff_polished.tsv | head -n1 > results/VNBI_genes.tsv
cut -f1,3,4,5,7,11,15,16,22 references/VNBI_liftoff.gff_polished.tsv | grep -e protein_coding_gene -e ncRNA_gene -e pseudogene >> results/VNBI_genes.tsv
 
head -n1 results/VNBI_genes.tsv > results/VNBI_loci_interest.tsv
#sed -i "s/$/\tLoci/" results/loci_interest.tsv

#awk '/SXI/,/STE12/' results/VNBI_genes.tsv >> results/VNBI_loci_interest.tsv
#sed -i "s/$/\tMAT/" results/loci_interest.tsv

cut -d',' -f2 results/centromeres.csv | while read line
do
awk -v i="$line" '$0 ~ i {print $0}' results/VNBI_genes.tsv
done >> results/VNBI_loci_interest.tsv