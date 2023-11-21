cd ./reference_genomes_test

#sed 's/protein_coding_gene/gene/g' FungiDB-65_CneoformansH99.gff > FungiDB-65_CneoformansH99_CZM.gff

grep -v "#" FungiDB-65_CneoformansH99.gff | cut -f 3 | sort | uniq > features.txt

#### FDB-65 vs. FDB-65 ####

liftoff -g FungiDB-65_CneoformansH99.gff -polish -o FungiDB-65to65_liftoff.gff -p 4 -f features.txt -u FDB65-FDB65_unmapped_features.txt FungiDB-65_CneoformansH99_Genome.fasta FungiDB-65_CneoformansH99_Genome.fasta

#### FDB-53 vs. FDB-53 ####

liftoff -g FungiDB-53_CneoformansH99.gff -polish -o FungiDB-53to53_liftoff.gff -p 4 -f features.txt -u FDB53-FDB53_unmapped_features.txt FungiDB-53_CneoformansH99_Genome.fasta FungiDB-53_CneoformansH99_Genome.fasta

#### FDB-53 vs. FDB-65 ####

liftoff -g FungiDB-53_CneoformansH99.gff -polish -o FungiDB-53to65_liftoff.gff -p 4 -f features.txt -u FDB53-FDB65_unmapped_features.txt FungiDB-65_CneoformansH99_Genome.fasta FungiDB-53_CneoformansH99_Genome.fasta

#### FDB-65 vs. FDB-53 ####

liftoff -g FungiDB-65_CneoformansH99.gff -polish -o FungiDB-65to53_liftoff.gff -p 4 -f features.txt -u FDB65-FDB53_unmapped_features.txt FungiDB-53_CneoformansH99_Genome.fasta FungiDB-65_CneoformansH99_Genome.fasta

# Count the number of appearances of each feature in every gff file

ls *gff* | grep -v "db"  | while read file
    do 
        cat features.txt | while read feature
            do 
                count=$(cut -f3 $file | grep ^$feature | wc -l)
                echo $file,$feature,$count
            done
    done > feature_count.csv