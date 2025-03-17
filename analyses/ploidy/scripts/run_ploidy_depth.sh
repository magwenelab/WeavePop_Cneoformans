
# Run ploidy depth analysis for all Ashton samples with putative duplications

grep "Desjardins" /FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/duplications_putative.tsv |\
cut -f8 |\
while read line
do
    quarto render ploidy_depth.qmd --output ploidy_depth_${line}.html -P sample_id:\"${line}\" -P dataset:"Desjardins" -P max_iterations:100000 -P threshold_weight:0.2 -P threshold_norm_depth:5
done
