
# Run the ploidy depth analysis for each sample in the metadata file 

mkdir -p analyses/notebooks/scripts/ploidy_depth/logs

for set in "Desjardins" "Ashton"
do
    grep $set analyses/results/tables/duplications_putative.tsv |\
    cut -f 8 |\
    while read line
    do
        quarto render analyses/scripts/ploidy_depth.qmd  -P sample_id:\"${line}\" -P dataset:\"${set}\" -P max_iterations:10 -P threshold_weight:0.2 -P threshold_norm_depth:5
        mv analyses/notebooks/scripts/ploidy_depth.html analyses/notebooks/scripts/ploidy_depth/${line}.html
    done
done
