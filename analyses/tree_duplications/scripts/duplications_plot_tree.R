# Use this script with the quarto environment
setwd("/FastData/czirion/Crypto_Diversity_Pipeline/")

library(tidyverse)
library(ggbeeswarm)
library(ggtree)
library(ggtreeExtra)
library(ape)
# library(phytools)
library(ggnewscale)
library(RColorBrewer)
#### Data
# Load the necessary data
metadata <- read.csv("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/data/derived/metadata_fixed.csv")
duplications <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/tables/polished/full_duplicated.tsv", sep = "\t", header = TRUE, stringsAsFactors = TRUE)

# Prepare the duplications_full data frame
duplications_full <- duplications %>%
    select(strain, chromosome) %>%
    distinct()

#### Plot the tree with duplicated chromosomes ####
# Make matrix of duplicated chromosomes
dup_chroms <- duplications_full %>%
    select(strain, chromosome)%>%
    mutate(duplicated_full = 1)%>%
    arrange(chromosome)%>%
    pivot_wider(names_from = chromosome, values_from = duplicated_full, values_fill = 0)%>%
    column_to_rownames("strain")%>%
    mutate(across(everything(), ~ ifelse(. == 1, cur_column(),"Euploid")))

euploid_strain <- metadata %>%
    filter(!strain %in% duplications_full$strain)%>%
    select(strain)

for (chrom in colnames(dup_chroms)){
    euploid_strain[chrom] <- "Euploid"
}

dup_chroms <- euploid_strain %>%
    column_to_rownames("strain") %>%
    bind_rows(dup_chroms)

# Get metadata

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
desj_tree_path <- "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/data/processed/desj_tree.newick"
desj_tree <- read.tree(desj_tree_path)

### Ashton tree
ashton_tree_path <- "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/data/processed/ashton_tree.newick"
ashton_tree <- read.tree(ashton_tree_path)

# Merged tree
merged_tree_path <- "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/data/processed/merged_tree.newick"
tree <- read.tree(merged_tree_path)

# Remove tips that are not in metadata$strain
tree <- drop.tip(tree, setdiff(tree$tip.label, metadata$strain))

# display colorblind friendly palettes
# display.brewer.all(colorblindFriendly = TRUE)
# display.brewer.all(colorblindFriendly = FALSE)

dataset_colors <- c("white", brewer.pal(9, "Set1")[c(1, 2)])
names(dataset_colors) <- levels(as.factor(metadata$dataset))

lineage_colors <- brewer.pal(8, "Dark2")[c(1, 2, 3, 4)]
names(lineage_colors) <- levels(as.factor(metadata$lineage))

sublineage_colors <- c("white", brewer.pal(12, "Set3")[c(1:10)])
names(sublineage_colors) <- levels(as.factor(metadata$vni_subdivision))

chrom_colors <- c(brewer.pal(nlevels(duplications$chromosome), "Paired"), "grey93")
names(chrom_colors) <- c(levels(duplications$chromosome), "Euploid")

source_colors <- brewer.pal(11, "BrBG")[c(9, 3)] # 9, 3 are the colors for the two sources
names(source_colors) <- levels(as.factor(metadata$source))



p <- ggtree(tree, layout = "circular", size = 0.1, branch.length = "none") + 
    geom_tiplab(aes(label = label), size = 0.2, align =TRUE, 
                    linetype = "dashed", linesize = .05)

p1 <- gheatmap(p, dataset, width=.05, colnames=FALSE, offset=3) +
    scale_fill_manual(values = dataset_colors, name="Dataset", na.translate = FALSE)+
    guides(fill = guide_legend(order = 1))+
    new_scale_fill()

p2 <- gheatmap(p1, lineage, width=.05, colnames=FALSE, offset=5) +
    scale_fill_manual(values = lineage_colors, name="Lineage", na.translate = FALSE)+
    guides(fill = guide_legend(order = 2))+
    new_scale_fill()

p3 <- gheatmap(p2, sublineage, width=.05, colnames=FALSE, offset=7,) +
    scale_fill_manual(values = sublineage_colors, name="VNI Sublineage", na.translate = FALSE)+
    guides(fill = guide_legend(order = 3))+
    new_scale_fill()

p4 <- gheatmap(p3, source, width=.05, colnames=FALSE, offset=9,) +
        scale_fill_manual(values = source_colors, name="Source", na.translate = FALSE)+
        guides(fill = guide_legend(order = 4))+
        new_scale_fill()

p5 <- gheatmap(p4, dup_chroms, width=.32, colnames = FALSE, offset=11,) +
    scale_fill_manual(values = chrom_colors, name="Duplicated\nchromosomes", na.translate = FALSE )+
    guides(fill = guide_legend(order = 5))+
    theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_text(size=9),
        legend.text=element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"))
p5
ggsave("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/figures/tree_merged_duplications.png", p5, height = 7, width = 7, units = "in", dpi = 900)

#### Plot tree with duplications of chromosomes 12 and 13 ####
dup_chroms_12_13 <- dup_chroms %>%
    select(chr12, chr13)

p5 <- gheatmap(p4, dup_chroms_12_13, width=.1, colnames = FALSE, offset=11,) +
    scale_fill_manual(values = chrom_colors, name="Duplicated\nchromosomes", na.translate = FALSE )+
    guides(fill = guide_legend(order = 5))+
    theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_text(size=9),
        legend.text=element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"))
p5
ggsave("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/figures/tree_merged_duplications_12_13.png", p5, height = 7, width = 7, units = "in", dpi = 900)

#### Plot the tree with only the samples that have duplications and the references####
keep_strains <- c(levels(duplications_full$strain), "H99", "Bt22", "Bt81")
tree_dups <- drop.tip(tree, setdiff(tree$tip.label, keep_strains))

p <- ggtree(tree_dups, layout = "rectangular", size = 0.5, branch.length = "none") + 
    geom_tiplab(aes(label = label), size = 3, align =TRUE, 
                    linetype = "dashed", linesize = 0.1, offset = 1)

p1 <- gheatmap(p, dataset, width=0.1, colnames=FALSE, offset=8) +
    scale_fill_manual(values = dataset_colors, name="Dataset", na.translate = FALSE)+
    guides(fill = guide_legend(order = 1))+
    new_scale_fill()

p2 <- gheatmap(p1, lineage, width=0.1, colnames=FALSE, offset=9.5) +
    scale_fill_manual(values = lineage_colors, name="Lineage", na.translate = FALSE)+
    guides(fill = guide_legend(order = 2))+
    new_scale_fill()

p3 <- gheatmap(p2, sublineage, width=0.1, colnames=FALSE, offset=11) +
    scale_fill_manual(values = sublineage_colors, name="VNI Sublineage", na.translate = FALSE)+ 
    guides(fill = guide_legend(order = 3))+
    new_scale_fill()

p4 <- gheatmap(p3, source, width=0.1, colnames=FALSE, offset=12.5) +
        scale_fill_manual(values = source_colors, name="Source", na.translate = FALSE)+
        guides(fill = guide_legend(order = 4))+
        new_scale_fill()

p5 <- gheatmap(p4, dup_chroms, width=0.7, colnames = FALSE, offset=14) +
    scale_fill_manual(values = chrom_colors, name="Duplicated\nchromosomes", na.translate = FALSE )+
    guides(fill = guide_legend(order = 5))+
    theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_text(size=9),
        legend.text=element_text(size=7),
        legend.key.size = unit(0.5, "cm"))

ggsave("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/figures/tree_merged_duplications_only_duplicated.png", p5, height = 7, width = 9, units = "in", dpi = 900)

p <- ggtree(tree_dups, layout = "rectangular", size = 0.5, branch.length = "none") + 
    geom_tiplab(aes(label = label), size = 3, align =TRUE, 
                    linetype = "dashed", linesize = 0.1, offset = 1)

p1 <- gheatmap(p, lineage, width=0.1, colnames=FALSE, offset=8) +
    scale_fill_manual(values = lineage_colors, name="Lineage", na.translate = FALSE)+
    guides(fill = guide_legend(order = 2))+
    new_scale_fill()

p2 <- gheatmap(p1, sublineage, width=0.1, colnames=FALSE, offset=10) +
    scale_fill_manual(values = sublineage_colors, name="VNI Sublineage", na.translate = FALSE)+ 
    guides(fill = guide_legend(order = 3))+
    new_scale_fill()

p3 <- gheatmap(p2, dup_chroms, width=0.7, colnames = FALSE, offset=12) +
    scale_fill_manual(values = chrom_colors, name="Duplicated\nchromosomes", na.translate = FALSE )+
    guides(fill = guide_legend(order = 5))+
    theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_text(size=9),
        legend.text=element_text(size=7),
        legend.key.size = unit(0.5, "cm"))

ggsave("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/figures/tree_merged_duplications_only_duplicated2.png", p3, height = 7, width = 9, units = "in", dpi = 900)

t<- ggtree(tree_dups, layout = "rectangular", size = 0.5, branch.length = "none") + 
    geom_tiplab(aes(label = label)) + 
    geom_nodelab(aes(label=node), geom="label")

d <- data.frame(node=c(104, 60, 57, 63), type=c("VNII", "VNBI", "VNBII", "VNI"))

p <- ggtree(tree_dups, layout = "rectangular", size = 0.5, branch.length = "none") + 
    geom_tiplab(aes(label = label), size = 3, align =TRUE, 
                    linetype = "dashed", linesize = 0.1, offset = 1)+
    geom_hilight(data=d, aes(node=node, fill=type), type = "rect")+
    scale_fill_manual(values = lineage_colors, name="Lineage", na.translate = FALSE)+
    guides(fill = guide_legend(order = 1))+
    new_scale_fill()

p2 <- gheatmap(p, dup_chroms, width=0.7, colnames = FALSE, offset=10) +
    scale_fill_manual(values = chrom_colors, name="Duplicated\nchromosomes", na.translate = FALSE )+
    guides(fill = guide_legend(order = 2))+
    new_scale_fill()

p3 <- gheatmap(p2, sublineage, width=0.1, colnames=FALSE, offset=19.5) +
    scale_fill_manual(values = sublineage_colors, name="VNI Sublineage", na.translate = FALSE)+ 
    guides(fill = guide_legend(order = 3))+
    theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_text(size=9),
        legend.text=element_text(size=7),
        legend.key.size = unit(0.5, "cm"))

ggsave("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/figures/tree_merged_duplications_only_duplicated3.png", p3, height = 7, width = 9, units = "in", dpi = 900)
