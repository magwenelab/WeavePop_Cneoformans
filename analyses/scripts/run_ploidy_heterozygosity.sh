# Run the ploidy heterozygosity analysis for each sample in the metadata file 
# (maybe I should change the file from which the samples are taken, so we don't analyze all of them)

for set in "Desjardins" "Ashton"
do
    grep $set analyses/data/processed/metadata_ashton_desj_all_weavepop_H99.csv |\
    cut -d',' -f1 |\
    while read line
    do
        quarto render analyses/scripts/ploidy_heterozygosity.qmd -P sample_id:\"${line}\" -P dataset:\"${set}\" 
        mkdir -p analyses/notebooks/scripts/ploidy_heterozygosity/
        mv analyses/notebooks/scripts/ploidy_heterozygosity.html analyses/notebooks/scripts/ploidy_heterozygosity/${line}.html
    done
done
