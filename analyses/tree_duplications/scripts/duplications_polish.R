# Use this script with the quarto environment
setwd("/FastData/czirion/Crypto_Diversity_Pipeline/")

library(tidyverse)

# After manual inspection of plots (see script duplications_gather_plots.xsh)
# I found that the following samples have partial duplications:
partial <- data.frame(sample= c("ERS2541051","ERS2541051","ERS542397","ERS1142798","ERS542490","ERS1142878","ERS2541358","SRS881238","ERS542397","ERS542498"),
                        chromosome = c("chr09","chr12","chr02","chr14","chr04","chr09","chr02","chr09","chr04","chr04"))

partial$sample_chromosome <- paste(partial$sample, partial$chromosome, sep="_")

# Remove partial duplications
duplications <- read_tsv("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/putative/duplications.tsv")

duplications_full <- duplications %>%
    mutate(sample_chromosome = paste(sample, chromosome, sep="_"))%>%
    filter(!sample_chromosome %in% partial$sample_chromosome)

duplicated_full_strain <- duplications_full %>%
    select(dataset, lineage, sample, strain, source, chromosome, norm_chrom_median, percent_cnv_size)

write_tsv(duplicated_full_strain, "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/polished/full_duplicated.tsv")

# Get multiple summary tables
dup_sample <- duplications_full%>%
    group_by(dataset,lineage, sample, strain, source) %>%
    summarise(n_chroms = n_distinct(chromosome), chromosomes = paste(chromosome, collapse = ", ")) %>%
    arrange(desc(n_chroms))
write_tsv(dup_sample, "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/polished/full_dup_sample.tsv")

dup_dataset_lineage_chromosome <-duplications_full%>%
    group_by(dataset,lineage, chromosome) %>%
    summarise(n_samples = n_distinct(sample), samples_in_dataset_lineage = first(samples_in_dataset_lineage))%>%
    mutate(percent_samples = round((n_samples / samples_in_dataset_lineage) * 100, 1))%>%
    select(dataset,lineage, chromosome, n_samples, samples_in_dataset_lineage, percent_samples)%>%
    arrange(chromosome, desc(lineage), desc(n_samples))

write_tsv(dup_dataset_lineage_chromosome, "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/polished/full_dup_dataset_lineage_chromosome.tsv")

dup_lineage_chromosome <- duplications_full %>%
    group_by(lineage, chromosome) %>%
    summarise(n_samples = n_distinct(sample), samples_in_lineage = first(samples_in_lineage))%>%
    mutate(percent_samples = round((n_samples / samples_in_lineage) * 100, 1))%>%
    select(lineage, chromosome, n_samples, samples_in_lineage,percent_samples)%>%
    arrange(chromosome, desc(lineage), desc(n_samples))

write_tsv(dup_lineage_chromosome, "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/polished/full_dup_lineage_chromosome.tsv")

dup_lineage_dataset <- duplications_full %>%
    group_by(dataset,lineage) %>%
    summarise(n_samples = n_distinct(sample), samples_in_dataset_lineage = first(samples_in_dataset_lineage))%>%
    mutate(percent_samples = round((n_samples / samples_in_dataset_lineage) * 100, 1))%>%
    select(lineage, n_samples, samples_in_dataset_lineage, percent_samples)%>%
    arrange(desc(lineage), desc(n_samples))

write_tsv(dup_lineage_dataset, "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/polished/full_dup_lineage_dataset.tsv")

dup_chromosome <- duplications_full %>%
    group_by(chromosome) %>%
    summarise(n_samples = n_distinct(sample), total_samples = first(total_samples))%>%
    mutate(percent_samples = round((n_samples / total_samples) * 100, 1))%>%
    select(chromosome, n_samples,total_samples, percent_samples)%>%
    arrange(chromosome, desc(n_samples))

write_tsv(dup_chromosome, "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/polished/full_dup_chromosome.tsv")

