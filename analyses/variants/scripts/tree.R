# install.packages("BiocManager")
# BiocManager::install("ggtree")
# BiocManager::install("ggtreeExtra")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ggtreeExtra))
suppressPackageStartupMessages(library(ape))
setwd("/FastData/czirion/Crypto_Desjardins/fungal_pop")

# =============================================================================
METADATA_FILE <- "/FastData/czirion/Crypto_Desjardins/config/metadata.csv"
metadata <- read.delim(METADATA_FILE, sep = ",", header = TRUE)
metadata <- metadata %>%
    select(strain, lineage)

strains_per_lineage <- metadata %>%
    group_by(lineage) %>%
    summarise(total_strains = n_distinct(strain))

metadata <- metadata %>%
    add_row(strain = "H99", lineage = "VNI")

TREE_FILE <- "CryptoDiversity_Desjardins_Tree.tre" 
# Note: The tree has H99 but the variants results don't
tree_des <- read.tree(TREE_FILE)

outgroup_clade <- metadata %>%
    filter(lineage == "VNII")%>%
    select(strain)

outgroup_clade <- outgroup_clade$strain
tree_des <- root(tree_des, outgroup = outgroup_clade)


# =============================================================================
#### ALL impact variants ####
ALL <- read_csv("ALL_filtered.csv", col_names = TRUE)
ALL$impact <- factor(ALL$impact, levels = c("HIGH","MODERATE", "LOW"))

# Filter out variants that are present in all strains of a lineage

strains_per_var <- ALL %>%
    group_by(lineage, var_id, impact) %>%
    summarise(num_strains = n_distinct(strain))

ALL_joined <- ALL %>%
    select(lineage, strain, var_id, gene_id, impact)%>%
    left_join(strains_per_var, by = c("var_id", "lineage", "impact"))%>%
    left_join(strains_per_lineage, by = "lineage")

ALL_filtered <- ALL_joined %>%
    filter(num_strains != total_strains)

# Number of genes with variants 

number_genes_vars <- ALL_filtered %>%
    group_by(strain, impact, gene_id) %>%
    summarise(num_vars = n_distinct(var_id)) %>%
    ungroup() %>%
    group_by(strain, impact) %>%
    summarise(num_genes = n_distinct(gene_id), mean_vars = mean(num_vars))

expanded <- expand.grid(
    strain = unique(metadata$strain),
    impact = unique(ALL_filtered$impact))

genes_per_strain <- left_join(expanded, number_genes_vars, by = c("strain", "impact"))
genes_per_strain$num_genes[is.na(genes_per_strain$num_genes)] <- 0

genes_per_strain <- left_join(genes_per_strain, metadata, by = "strain")

# =============================================================================
# Add genes with variants to tree


genes_per_strain_wide <- genes_per_strain %>%
    select(strain, num_genes, impact)%>%
    pivot_wider(names_from = impact, values_from = num_genes, values_fill = 0)

plot_tree <- ggtree(tree_des) %<+% metadata + 
                geom_tiplab(aes(label = label), 
                size = 2, align = TRUE)

plot_tree
ggsave("plot_tree.png", plot_tree, height = 20, width = 10, units = "in", dpi = 300)
var_tree <- plot_tree + 
            geom_facet(panel= "HIGH",
                        data = genes_per_strain_wide,
                        geom = geom_col,
                        mapping = aes(x = HIGH,color = lineage, fill = lineage),
                        orientation = "y") +
            geom_facet(panel= "MODERATE",
                        data = genes_per_strain_wide,
                        geom = geom_col,
                        mapping = aes(x = MODERATE,fill =lineage, color = lineage),
                        orientation = "y") +
            geom_facet(panel= "LOW",
                        data = genes_per_strain_wide,
                        geom = geom_col,
                        mapping = aes(x = LOW,fill = lineage, color = lineage),
                        orientation = "y") +
            scale_y_continuous(expand=c(0, 0.3)) +
            theme_tree2(legend.position=c(.05, .85))

var_tree
ggsave("var_tree.png", var_tree, height = 20, width = 10, units = "in", dpi = 300)

