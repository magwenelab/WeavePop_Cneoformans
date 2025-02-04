# Use this script with the quarto environment
setwd("/FastData/czirion/Crypto_Diversity_Pipeline/")

library(tidyverse)
library(ggbeeswarm)
library(ggtree)
library(ggtreeExtra)
library(ape)
library(phytools)
library(ggnewscale)
library(RColorBrewer)

#### Metadata #### ========================================================================================================================
metadata <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/data/derived/metadata_fixed.csv", header=TRUE, sep=",")

metadata <- metadata %>%
    select(sample, strain, source, lineage, dataset,name, vni_subdivision)%>%
    filter(!strain == "H99") %>%
    mutate(source = ifelse(source == "Environmental", "Environment", source))%>%
    group_by(dataset, lineage)%>%
    mutate(samples_in_dataset_lineage = n_distinct(sample))%>%
    ungroup() %>%
    group_by(lineage)%>%
    mutate(samples_in_lineage = n_distinct(sample))%>%
    ungroup()%>%
    mutate(total_samples = n_distinct(sample))

chromosome_names = read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Desjardins_Ashton/results_joined_241204/02.Dataset/chromosomes.csv", header=TRUE, sep=",")
chromosome_names <- chromosome_names %>%
    mutate(chromosome = str_pad(chromosome, 2, pad = "0"))%>%
    mutate(chromosome = as.factor(chromosome))
levels(chromosome_names$chromosome) <- paste("chr", chromosome_names$chromosome, sep="")

#### Get metrics of chromosomes from called CNVs #### ==========================================================================

cnv_calls <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Desjardins_Ashton/results_joined_241204/02.Dataset/cnv/cnv_calls.tsv", header=TRUE , sep="\t")
chromosome_lengths = read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/data/derived/chromosome_lengths.tsv", header=FALSE, col.names=c("accession", "length"), sep="\t")

repeat_fraction_threshold <- 0.2

# Percentage of each chromosome duplicated 
cnv_percent <- cnv_calls %>%
    filter(repeat_fraction < repeat_fraction_threshold) %>%
    filter(cnv == "duplication")%>%
    group_by(accession, sample, cnv) %>%
    summarise(total_cnv_size = sum(region_size),
                n_cnvs = n(),
                first = min(start),
                last = max(end),
                mean_smooth_depth = round(mean(smooth_depth),2)) %>%
    left_join(chromosome_lengths, by="accession") %>%
    left_join(chromosome_names, by="accession") %>%
    left_join(metadata, by=c("sample", "lineage")) %>%
    mutate(percent_cnv_size = round((total_cnv_size / length) * 100, 2),
            size_covered = last - first,
            percent_size_covered = round((size_covered / length) * 100, 2))%>%
    mutate(chromosome = as.factor(chromosome)) %>%
    select(dataset,lineage, samples_in_lineage, samples_in_dataset_lineage,total_samples, sample,strain, source, accession, chromosome, percent_cnv_size, percent_size_covered, mean_smooth_depth, n_cnvs)


boxplot <- ggplot(cnv_percent)+
    geom_quasirandom(aes(x=chromosome, y=percent_cnv_size, color = source))+
    facet_grid(lineage~chromosome, scales="free")+
    labs(x="Chromosome", y="Percent of chromosome duplicated (%)", title=paste("Distribution of the size of duplications with ", repeat_fraction_threshold, "fraction of repeats allowed"))+
    ylim(0, 100)+
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "lightgray"),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2))
ggsave("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/figures/boxplot_cnv_percent.png", boxplot, height = 8, width = 8, units = "in", dpi = 500)

#### Percentage of each chromosome covered by repeats #### ==========================================================================

lineages <- unique(metadata$lineage)
repeats_all <- data.frame()
for(lineage in lineages){
    repeats_path <-paste("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Desjardins/results_241202/03.References/", lineage, "/", lineage, "_repeats.bed", sep ="")
    repeats <- read.delim(repeats_path, header=FALSE, col.names=c("accession", "start", "end", "repeat_type"), sep="\t")
    repeats$lineage <- lineage
    repeats_all <- rbind(repeats_all, repeats)
}

repeats_percent <- repeats_all %>%
    mutate(repeat_size_each = end - start)%>%
    group_by(accession, lineage) %>%
    summarise(repeat_size = sum(repeat_size_each)) %>%
    left_join(chromosome_lengths, by="accession") %>%
    left_join(chromosome_names, by=c("accession","lineage")) %>%
    mutate(percent_repeat_size = round((repeat_size / length) * 100, 2))%>%
    mutate(chromosome = as.factor(chromosome))%>%
    select(lineage, accession, chromosome, percent_repeat_size)

barplot <- ggplot(repeats_percent)+
    geom_col(aes(x=lineage, y=percent_repeat_size, fill = lineage))+
    facet_grid(~chromosome, scales="free")+
    labs(title = "Percentage of each chromosome covered by repeats", x="Lineage", y="Percent of chromosome covered by repeats (%)")+
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "lightgray"),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2))
ggsave("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/figures/barplot_repeats.png", barplot, height = 8, width = 8, units = "in", dpi = 500)   

#### Get chromosome median depth #### ==========================================================================
depth_by_chrom_good_desjardins <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Desjardins/results_241202/04.Intermediate_files/02.Dataset/depth_quality/depth_by_chrom_good.tsv", header=TRUE, sep="\t")
depth_by_chrom_good_ashton <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Ashton/results_241203/04.Intermediate_files/02.Dataset/depth_quality/depth_by_chrom_good.tsv", header=TRUE, sep="\t")
depth_by_chrom_good <- rbind(depth_by_chrom_good_desjardins, depth_by_chrom_good_ashton)
depth_by_chrom <- depth_by_chrom_good %>%
    select(sample, accession, norm_chrom_median)


#### Join CNV metrics with normalized chromosome depth #### ==========================================================================

cnv_and_depth <- left_join(depth_by_chrom, cnv_percent, by = c("sample", "accession" ))

percent_size_threshold <- 80
depth_threshold <- 1.55

# Filter chromosomes by normalized depth OR percent of the chromosome in CNV regions
duplications <- cnv_and_depth %>%
    filter(norm_chrom_median > depth_threshold | percent_cnv_size > percent_size_threshold)%>%
    select(dataset, lineage, samples_in_lineage, samples_in_dataset_lineage, total_samples, source, strain, sample, chromosome, norm_chrom_median, mean_smooth_depth, percent_cnv_size, percent_size_covered, n_cnvs)


# Plot all metrics of putative duplications
duplications_long <- duplications  %>%
    arrange(percent_cnv_size) %>%
    mutate(chrom_sample = paste(chromosome, sample, sep="_"))
duplications_long$chrom_sample <- factor(duplications_long$chrom_sample, levels = duplications_long$chrom_sample)
duplications_long <- pivot_longer(duplications_long, cols = c("percent_cnv_size", "percent_size_covered", "norm_chrom_median", "mean_smooth_depth", "n_cnvs"), names_to = "variable", values_to = "value")

p <- ggplot(duplications_long, aes(x = chrom_sample, y = value, color = variable))+
        geom_point()+
        facet_grid(scales="free", rows= vars(variable))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

#### Get multiple summary tables ==========================================================================

write_tsv(duplications, "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/putative/duplications.tsv")


dup_sample <- duplications %>%
    group_by(dataset,lineage, sample, strain, source) %>%
    summarise(n_chroms = n_distinct(chromosome), chromosomes = paste(chromosome, collapse = ", ")) %>%
    arrange(desc(n_chroms))
write_tsv(dup_sample, "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/putative/dup_sample.tsv")

dup_dataset_lineage_chromosome <- duplications %>%
    group_by(dataset,lineage, chromosome) %>%
    summarise(n_samples = n_distinct(sample), samples_in_dataset_lineage = first(samples_in_dataset_lineage))%>%
    mutate(percent_samples = round((n_samples / samples_in_dataset_lineage) * 100, 1))%>%
    select(dataset,lineage, chromosome, n_samples, samples_in_dataset_lineage, percent_samples)%>%
    arrange(chromosome, desc(lineage), desc(n_samples))

write_tsv(dup_dataset_lineage_chromosome, "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/putative/dup_dataset_lineage_chromosome.tsv")

dup_lineage_chromosome <- duplications%>%
    group_by(lineage, chromosome) %>%
    summarise(n_samples = n_distinct(sample), samples_in_lineage = first(samples_in_lineage))%>%
    mutate(percent_samples = round((n_samples / samples_in_lineage) * 100, 1))%>%
    select(lineage, chromosome, n_samples, samples_in_lineage,percent_samples)%>%
    arrange(chromosome, desc(lineage), desc(n_samples))

write_tsv(dup_lineage_chromosome, "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/putative/dup_lineage_chromosome.tsv")

dup_lineage_dataset <- duplications%>%
    group_by(dataset,lineage) %>%
    summarise(n_samples = n_distinct(sample), samples_in_dataset_lineage = first(samples_in_dataset_lineage))%>%
    mutate(percent_samples = round((n_samples / samples_in_dataset_lineage) * 100, 1))%>%
    select(lineage, n_samples, samples_in_dataset_lineage, percent_samples)%>%
    arrange(desc(lineage), desc(n_samples))

write_tsv(dup_lineage_dataset, "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/putative/dup_lineage_dataset.tsv")

dup_chromosome <- duplications %>%
    group_by(chromosome) %>%
    summarise(n_samples = n_distinct(sample), total_samples = first(total_samples))%>%
    mutate(percent_samples = round((n_samples / total_samples) * 100, 1))%>%
    select(chromosome, n_samples,total_samples, percent_samples)%>%
    arrange(chromosome, desc(n_samples))

write_tsv(dup_chromosome, "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/putative/dup_chromosome.tsv")

