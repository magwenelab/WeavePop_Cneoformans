# Use this script with the quarto environment
setwd("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/scripts")

library(tidyverse)
library(ggbeeswarm)
library(ggtree)
library(ggtreeExtra)
library(ape)
library(phytools)
library(ggnewscale)

#### Metadata ####
metadata <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/data/metadata_fixed.csv", header=TRUE, sep=",")

#### Duplicated samples from chromosome median depth ####
depth_by_chrom_good_desjardins <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Desjardins/results_241202/04.Intermediate_files/02.Dataset/depth_quality/depth_by_chrom_good.tsv", header=TRUE, sep="\t")
depth_by_chrom_good_ashton <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Ashton/results_241203/04.Intermediate_files/02.Dataset/depth_quality/depth_by_chrom_good.tsv", header=TRUE, sep="\t")
depth_by_chrom_good <- rbind(depth_by_chrom_good_desjardins, depth_by_chrom_good_ashton)
chromosome_names = read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/results_joined_241204/02.Dataset/chromosomes.csv", header=TRUE, sep=",")

depth_threshold <- 1.55

metadata <- metadata %>%
    select(sample, strain, source, lineage, dataset,name, vni_subdivision)%>%
    group_by(dataset, lineage)%>%
    mutate(samples_per_lineage = n_distinct(sample))%>%
    ungroup() %>%
    mutate(source = ifelse(source == "Environmental", "Environment", source))

chromosome_names <- chromosome_names %>%
    mutate(chromosome = str_pad(chromosome, 2, pad = "0"))%>%
    mutate(chromosome = as.factor(chromosome))
levels(chromosome_names$chromosome) <- paste("chr", chromosome_names$chromosome, sep="")

# Filter samples with median depth above threshold
duplicated <- depth_by_chrom_good %>%
    filter(norm_chrom_median > depth_threshold)%>%
    left_join(metadata, by="sample")%>%
    left_join(chromosome_names, by=c("accession", "lineage"))%>%
    select(dataset,lineage,samples_per_lineage,sample,strain, source, accession, chromosome, norm_chrom_median)

write_tsv(duplicated, "../results/duplicated_samples.tsv")

multiple_duplication <- duplicated %>%
    group_by(dataset,lineage, sample, strain, source) %>%
    summarise(n_chroms = n_distinct(chromosome), chromosomes = paste(chromosome, collapse = ", "))

duplicated_summary <- duplicated %>%
    group_by(dataset,lineage, chromosome) %>%
    summarise(n_samples = n_distinct(sample), samples_per_lineage = first(samples_per_lineage))%>%
    mutate(percent_samples = round((n_samples / samples_per_lineage) * 100, 1))%>%
    select(dataset,lineage, chromosome, n_samples)%>%
    arrange(chromosome, desc(lineage), desc(n_samples))

write_tsv(duplicated_summary, "../results/duplicated_summary.tsv")

#### Plot the tree with duplicated chromosomes ####
dup_chroms <- duplicated %>%
    select(strain, chromosome)%>%
    mutate(duplicated = 1)%>%
    arrange(chromosome)%>%
    pivot_wider(names_from = chromosome, values_from = duplicated, values_fill = 0)%>%
    column_to_rownames("strain")%>%
    mutate(across(everything(), ~ ifelse(. == 1, cur_column(), NA)))

lineage <- metadata %>%
    select(strain, lineage)%>%
    column_to_rownames("strain")

source <- metadata %>%
    select(strain, source)%>%
    column_to_rownames("strain")

sublineage <- metadata %>%
    select(strain, vni_subdivision)%>%
    column_to_rownames("strain")

dataset <- metadata %>%
    select(strain, dataset)%>%
    column_to_rownames("strain")

# Desjardins tree
desj_tree_path <- "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/data/desjardins_tree.newick"
desj_tree <- read.tree(desj_tree_path)

### Ashton tree
ashton_tree_path <- "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/data/ashton_tree.newick"
ashton_tree <- read.tree(ashton_tree_path)

# Merged tree
merged_tree_path <- "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/data/merged_tree.newick"
tree <- read.tree(merged_tree_path)


p <- ggtree(tree, layout = "circular", size = 0.1) + 
    geom_tiplab(aes(label = label), size = 0.25, align =TRUE, 
                    linetype = "dashed", linesize = .05)+
    geom_treescale(x=0.6, y=0, width=0.01, offset = 4)s

p1 <- gheatmap(p, dataset, width=.05, colnames=FALSE, offset=.025) +
    scale_fill_brewer(palette = "Set1", name="Dataset",  na.translate = FALSE)+ 
    new_scale_fill()

p2 <- gheatmap(p1, lineage, width=.05, colnames=FALSE, offset=.042) +
    scale_fill_brewer(palette = "Dark2", name="Lineage",  na.translate = FALSE)+ 
    new_scale_fill()

p3 <- gheatmap(p2, sublineage, width=.05, colnames=FALSE, offset=.059,) +
    scale_fill_brewer(palette="Set3", name="VNI Sublineage", na.translate = FALSE)+ 
    new_scale_fill() 

p4 <- gheatmap(p3, source, width=.05, colnames=FALSE, offset=.076,) +
        scale_fill_brewer(palette="Set2", name="Source", na.translate = FALSE)+ 
        new_scale_fill()

p5 <- gheatmap(p4, dup_chroms, width=.32, colnames = FALSE, offset=0.095,) +
    scale_fill_brewer(palette="Paired", name="Duplicated\nchromosomes",
        na.value = "grey93")


ggsave("../results/duplications_merged_tree.png", p5, height = 7, width = 7, units = "in", dpi = 900)
ggsave("../results/duplications_merged_tree.svg", p5, height = 8, width = 8, units = "in")


#### Desjardins samples in Ashton ####
ashton_strains <- ashton_tree$tip.label
desjardins_in_ashton <- metadata %>%
    filter(strain %in% ashton_strains)%>%
    filter(dataset == "Desjardins")
nrow(desjardins_in_ashton)

desjardins_not_in_ashton <- metadata %>%
    filter(dataset == "Desjardins")%>%
    filter(lineage == "VNI")%>%
    filter(!strain %in% ashton_strains)

nrow(desjardins_not_in_ashton)

#### Get duplication from called CNVs ####

cnv <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/results_joined_241204/02.Dataset/cnv/cnv_calls.tsv", header=TRUE , sep="\t")
chromosome_lengths = read.delim("../data/chromosome_lengths.tsv", header=FALSE, col.names=c("accession", "length"), sep="\t")

repeat_fraction_threshold <- 0.2

#### Percentage of each chromosome duplicated ####
cnv_percent <- cnv %>%
    filter(repeat_fraction < repeat_fraction_threshold) %>%
    filter(cnv == "duplication")%>%
    group_by(accession, sample, cnv) %>%
    summarise(total_cnv_size = sum(region_size)) %>%
    left_join(chromosome_lengths, by="accession") %>%
    left_join(chromosome_names, by="accession") %>%
    left_join(metadata, by=c("sample", "lineage")) %>%
    mutate(percent_cnv_size = round((total_cnv_size / length) * 100, 1))%>%
    mutate(chromosome = as.factor(chromosome))

    
boxplot <- ggplot(cnv_percent)+
    geom_quasirandom(aes(x=chromosome, y=percent_cnv_size, color = source))+
    facet_grid(lineage~chromosome, scales="free")+
    labs(x="Chromosome", y="Percent of chromosome duplicated (%)", title=paste("Duplication size distribution with ", repeat_fraction_threshold, "fraction of repeats allowed"))+
    ylim(0, 100)+
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "lightgray"),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2))
ggsave("../results/cnv_boxplot.png", boxplot, height = 8, width = 8, units = "in", dpi = 500)

#### Percentage of each chromosome covered by repeats ####

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
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "lightgray"),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2))
ggsave("../results/repeats_barplot.png", barplot, height = 8, width = 8, units = "in", dpi = 500)   


#### Number of samples with duplicated chromosomes ####

percent_size_threshold <- 80

cnv_duplicated <- cnv_percent %>%
    filter(percent_cnv_size > percent_size_threshold)%>%
    select(dataset, lineage, samples_per_lineage, sample, strain, source, accession, chromosome, percent_cnv_size)

write_tsv(cnv_duplicated, "../results/cnv_duplicated_samples.tsv")
    
cnv_summary <- cnv_duplicated%>%
    group_by(dataset,lineage, chromosome) %>%
    summarise(n_samples = n_distinct(sample))%>%
    arrange(chromosome, desc(lineage), desc(n_samples))

write_tsv(cnv_summary, "../results/cnv_duplicated_summary.tsv")
