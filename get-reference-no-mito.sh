seqkit grep -n -r -v -p "mitochondrion" reference_genomes/VNBI.fasta > reference_genomes/VNBI.fasta.modif
mv reference_genomes/VNBI.fasta reference_genomes/VNBI.fasta.original
mv reference_genomes/VNBI.fasta.modif reference_genomes/VNBI.fasta