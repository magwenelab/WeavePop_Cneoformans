# Use this script with the quarto environment
setwd("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/scripts")

library(tidyverse)
library(ggbeeswarm)
library(ggtree)
library(ggtreeExtra)
library(ape)
library(phytools)
library(ggnewscale)


depth_by_chrom_good_desjardins <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Desjardins/results_241118/02.Dataset/depth_quality/depth_by_chrom_good.tsv", header=TRUE, sep="\t")
depth_by_chrom_good_ashton <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Ashton/results_241119/02.Dataset/depth_quality/depth_by_chrom_good.tsv", header=TRUE, sep="\t")
depth_by_chrom_good <- rbind(depth_by_chrom_good_desjardins, depth_by_chrom_good_ashton)

depth_threshold <- 1.55

metadata <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/results_joined_241120/02.Dataset/metadata.csv", header=TRUE, sep=",")
metadata <- metadata %>%
    select(sample, strain, source, lineage)%>%
    group_by(lineage)%>%
    mutate(samples_per_lineage = n_distinct(sample))%>%
    ungroup()

chromosome_names = read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/results_joined_241120/02.Dataset/chromosomes.csv", header=TRUE, sep=",")
chromosome_names <- chromosome_names %>%
    mutate(chromosome = str_pad(chromosome, 2, pad = "0"))%>%
    mutate(chromosome = as.factor(chromosome))
    
levels(chromosome_names$chromosome) <- paste("chr", chromosome_names$chromosome, sep="")

duplicated <- depth_by_chrom_good %>%
    filter(norm_chrom_median > depth_threshold)%>%
    left_join(metadata, by="sample")%>%
    left_join(chromosome_names, by=c("accession", "lineage"))%>%
    select(lineage,samples_per_lineage,sample,strain, source, accession, chromosome, norm_chrom_median)

multiple_duplication <- duplicated %>%
    group_by(lineage, sample, strain, source) %>%
    summarise(n_chroms = n_distinct(chromosome), chromosomes = paste(chromosome, collapse = ", "))

duplicated_summary <- duplicated %>%
    group_by(lineage, chromosome) %>%
    summarise(n_samples = n_distinct(sample), samples_per_lineage = first(samples_per_lineage))%>%
    arrange(desc(n_samples))%>%
    mutate(percent_samples = round((n_samples / samples_per_lineage) * 100, 1))%>%
    select(lineage, chromosome, n_samples, percent_samples)

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

tree_path <- "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/data/CryptoDiversity_Desjardins_Tree.tre"
tree <- read.tree(tree_path)

# outgroup_clade <- metadata %>%
#     filter(lineage == "VNII")%>%
#     select(strain)

# outgroup_clade <- outgroup_clade$strain
# tree <- root(tree, outgroup = outgroup_clade)
# tree <- root(tree, node = 493)
# plot(tree, show.tip.label = FALSE) 
# nodelabels()

# Reroot the tree at the middle of the branch leading to VNII
VNII_root <- getMRCA(tree, c("C2","C12"))

edge_length <- subset(tree$edge.length, tree$edge[,2] == VNII_root)
tree <- reroot(tree, VNII_root, edge_length/2)

p <- ggtree(tree, layout = "circular") + 
    geom_tiplab(aes(label = label), size = 0.8, align =TRUE, 
                    linetype = "dashed", linesize = .05)

p1 <- gheatmap(p, lineage, width=.08, colnames=FALSE, offset=.08) +
    scale_fill_brewer(palette = "Paired", name="Lineage",  na.translate = FALSE)

p2 <- p1 + new_scale_fill() 

p3 <- gheatmap(p2, source,  offset=.1, width=.08, colnames=FALSE) +
        scale_fill_brewer(palette="Set2", name="Source", na.translate = FALSE)

p4 <- p3 + new_scale_fill() 

p5 <- gheatmap(p4, dup_chroms,offset=0.12, width=.32, 
    colnames_offset_y = .25, colnames_angle = 90, font.size = 1.5) +
    scale_fill_brewer(palette="Dark2", name="Duplicated\nchromosomes",
        na.value = "grey93")
    # theme_tree(legend.position = "bottom", legend.direction = "vertical")
ggsave("cnv_tree.png", p5, height = 8, width = 8, units = "in", dpi = 500)

