setwd("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/scripts")

library(tidyverse)
library(ggbeeswarm)
library(ggtree)
library(ggtreeExtra)
library(ape)

cnv <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Ashton/results_joined_2024-10-16/02.Dataset/cnv/cnv_calls.tsv", header=TRUE , sep="\t")
chromosome_lengths = read.delim("../data/chromosome_lengths.tsv", header=FALSE, col.names=c("Accession", "length"), sep="\t")
chromosome_names = read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Desjardins/config/chromosomes.csv", header=TRUE, sep=",")
metadata <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Ashton/results_joined_2024-10-16/02.Dataset/metadata.csv", header=TRUE, sep=",")

metadata <- metadata %>%
    select("Sample" = sample, strain, source, "Lineage" = lineage)

metadata <- metadata %>%
    mutate(source = ifelse(source == "Environmental", "Environment", source))
colnames(chromosome_names) <- str_to_title(colnames(chromosome_names)) #

repeat_fraction_threshold <- 0.2

#### Percentage of each chromosome duplicated ####
cnv_percent <- cnv %>%
    filter(Repeat_fraction < repeat_fraction_threshold) %>%
    filter(CNV == "DUPLICATION")%>%
    group_by(Accession, Sample, CNV) %>%
    summarise(total_cnv_size = sum(Region_Size)) %>%
    left_join(chromosome_lengths, by="Accession") %>%
    left_join(chromosome_names, by="Accession") %>%
    left_join(metadata, by=c("Sample", "Lineage")) %>%
    mutate(percent_cnv_size = round((total_cnv_size / length) * 100, 1))%>%
    mutate(Chromosome = as.factor(Chromosome))

    
boxplot <- ggplot(cnv_percent)+
    geom_quasirandom(aes(x=Chromosome, y=percent_cnv_size, color = source))+
    facet_grid(Lineage~Chromosome, scales="free")+
    labs(x="Chromosome", y="Percent of chromosome duplicated (%)", title=paste("Duplication size distribution with ", repeat_fraction_threshold, "fraction of repeats allowed"))+
    ylim(0, 100)+
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "lightgray"),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2))
boxplot    

#### Percentage of each chromosome covered by repeats ####

lineages <- unique(metadata$Lineage)
repeats_all <- data.frame()
for(lineage in lineages){
    repeats_path <-paste("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Desjardins/results_241118/03.References/", lineage, "/", lineage, "_repeats.bed", sep ="")
    repeats <- read.delim(repeats_path, header=FALSE, col.names=c("Accession", "Start", "End", "Repeat_type"), sep="\t")
    repeats$lineage <- lineage
    repeats_all <- rbind(repeats_all, repeats)
}
repeats_percent <- repeats_all %>%
    mutate(repeat_size_each = End - Start)%>%
    group_by(Accession, lineage) %>%
    summarise(repeat_size = sum(repeat_size_each)) %>%
    left_join(chromosome_lengths, by="Accession") %>%
    left_join(chromosome_names, by="Accession") %>%
    mutate(percent_repeat_size = round((repeat_size / length) * 100, 2))%>%
    mutate(Chromosome = as.factor(Chromosome))%>%
    select(Lineage, Chromosome, percent_repeat_size)

barplot <- ggplot(repeats_percent)+
    geom_col(aes(x=Lineage, y=percent_repeat_size, fill = Lineage))+
    facet_grid(~Chromosome, scales="free")+
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "lightgray"),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2))
barplot    


#### Number of samples with duplicated chromosomes ####

percent_size_threshold <- 80
samples_per_lineage <- metadata %>%
    group_by(Lineage) %>%
    summarise(samp_per_lin = n_distinct(Sample))

cnv_summary <- cnv_percent %>%
    filter(percent_cnv_size > percent_size_threshold)%>%
    group_by(Lineage, Chromosome) %>%
    summarise(n_samples = n_distinct(Sample))%>%
    arrange(desc(n_samples))%>%
    left_join(samples_per_lineage, by="Lineage")%>%
    mutate(percent_samples = round((n_samples / samp_per_lin) * 100, 1))%>%
    select(Lineage, Chromosome, n_samples, percent_samples)


#### Using depth_by_chrom_good ####
depth_by_chrom_good <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Desjardins/results_241118/02.Dataset/depth_quality/depth_by_chrom_good.tsv", header=TRUE, sep="\t")
metadata <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Ashton/results_joined_2024-10-16/02.Dataset/metadata.csv", header=TRUE, sep=",")
metadata <- metadata %>%
    select(sample, strain, source, lineage)%>%
    group_by(lineage)%>%
    mutate(samples_per_lineage = n_distinct(sample))%>%
    ungroup()

chromosome_names = read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Desjardins/config/chromosomes.csv", header=TRUE, sep=",")
chromosome_names <- chromosome_names %>%
    mutate(chromosome = as.factor(chromosome))
    
levels(chromosome_names$chromosome) <- paste("chr", chromosome_names$chromosome, sep="")

duplicated <- depth_by_chrom_good %>%
    filter(norm_chrom_median > 1.55)%>%
    left_join(metadata, by="sample")%>%
    left_join(chromosome_names, by=c("accession", "lineage"))%>%
    select(lineage,samples_per_lineage,sample,strain, source, accession, chromosome, norm_chrom_median)

duplicated_summary <- duplicated %>%
    group_by(lineage, chromosome) %>%
    summarise(n_samples = n_distinct(sample), samples_per_lineage = first(samples_per_lineage))%>%
    arrange(desc(n_samples))%>%
    mutate(percent_samples = round((n_samples / samples_per_lineage) * 100, 1))%>%
    select(lineage, chromosome, n_samples, percent_samples)

df_info <- duplicated %>%
    select(strain, chromosome)%>%
    mutate(duplicated = 1)%>%
    arrange(chromosome)%>%
    pivot_wider(names_from = chromosome, values_from = duplicated, values_fill = 0)%>%
    column_to_rownames("strain")

df_info <- df_info %>%
    mutate(across(everything(), ~ ifelse(. == 1, cur_column(), NA)))




tree_path <- "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/data/CryptoDiversity_Desjardins_Tree.tre"
tree <- read.tree(tree_path)

outgroup_clade <- metadata %>%
    filter(lineage == "VNII")%>%
    select(strain)

outgroup_clade <- outgroup_clade$strain
tree <- root(tree, outgroup = outgroup_clade)


p <- ggplot(tree, aes(x, y)) + 
    geom_tree() + 
    geom_tiplab(aes(label = label), size = 2) +
    theme_tree()


p <- ggtree(tree, layout = "circular") + 
    # geom_tiplab(aes(label = label), size = 2) +
    theme_tree()

gheatmap(p, df_info, offset=.8, width=.2,
               colnames_angle=95, colnames_offset_y = .25) +
    scale_fill_viridis_d(option="D", name="Duplicated\nchromosome")

gheatmap(p, df_info,width=.2,) +
    scale_fill_viridis_d(option="D", name="Duplicated\nchromosome")
