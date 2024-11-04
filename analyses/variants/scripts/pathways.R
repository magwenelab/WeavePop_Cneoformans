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
suppressPackageStartupMessages(library(ggbeeswarm))
setwd("/FastData/czirion/Crypto_Desjardins/fungal_pop")

#### Pathway variants ####
create_barplot <- function(PWY) {
    PWY_name <- deparse(substitute(PWY))
    PWY$impact <- factor(PWY$impact, levels = c("HIGH","MODERATE", "LOW"))

    PWY_var_per_strain <- PWY %>%
        group_by(strain, gene_name,impact) %>%
        summarise(num_vars = n_distinct(var_id))

    i_colors <- brewer.pal(3, "Set1")

    plot <- ggplot(PWY_var_per_strain, aes(x = strain, y = num_vars, fill = impact )) +
        geom_col() +
        facet_grid(gene_name~lineage, scales = "free_x", space = "free_x") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
              legend.position = "top") +
        labs(title = paste("Number of variants per strain in the ", PWY_name, " pathway", sep = ""),
             x = "Strain",
             y = "Number of variants",
             fill = "Impact") +
        scale_fill_manual(values = i_colors)

    ggsave(paste("results/fungal_pop/", PWY_name, "_variants.png", sep = ""), plot, width = 18, height = 9)
}

create_boxplot<- function(PWY, metadata){
    PWY_name <- deparse(substitute(PWY))

    strains_per_lineage <- metadata %>%
        group_by(lineage) %>%
        summarise(total_strains = n_distinct(strain))

    strain_per_var <- PWY %>%
        group_by(var_id) %>%
        summarise(num_strains = n_distinct(strain))

    PWY_joined <- PWY %>%
        select(lineage, strain, var_id, impact, gene_name)%>%
        left_join(strain_per_var, by = "var_id")%>%
        left_join(strains_per_lineage, by = "lineage")%>%
        filter(num_strains != total_strains)

    PWY_joined$impact <- factor(PWY_joined$impact, levels = c("HIGH","MODERATE", "LOW"))

    number_vars <- PWY_joined %>%
        group_by(strain, gene_name, impact) %>%
        summarise(num_vars = n_distinct(var_id))
    
    expanded <- expand.grid(
        strain = unique(metadata$strain),
        gene_name = unique(PWY_joined$gene_name),
        impact = unique(PWY_joined$impact))

    var_per_strain <- left_join(expanded, number_vars, by = c("strain", "gene_name", "impact"))
    var_per_strain$num_vars[is.na(var_per_strain$num_vars)] <- 0
    
    var_per_strain <- left_join(var_per_strain, metadata, by = "strain")

    i_colors <- brewer.pal(3, "Set1")

    plot <- ggplot(var_per_strain, aes(x = lineage, y = num_vars)) +
        geom_quasirandom(aes(color = impact)) +
        geom_boxplot(alpha = 0) +
        facet_grid(impact~gene_name, scales = "free_y") +
        theme(panel.grid = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2),
            legend.position = "none") +
        labs(title = paste("Number of variants, separated by impact, per strain in the genes of the ", PWY_name, " pathway", sep = ""),
            x = "Lineage",
            y = "Number of variants") +
        scale_color_manual(values = i_colors) +
        scale_y_continuous(labels = scales::comma, limits = c(0, NA))
    ggsave(paste("results/fungal_pop/", PWY_name, "_boxplot.png", sep = ""), plot, width = 18, height = 9)
    
}

metadata <- read_csv("/FastData/czirion/Crypto_Desjardins/config/metadata.csv", col_names = TRUE)
metadata <- metadata %>%
    select(strain, lineage)

TOR <- read_csv("results/fungal_pop/TOR.csv", col_names = TRUE)
create_barplot(TOR)
create_boxplot(TOR, metadata)

HOG <- read_csv("results/fungal_pop/HOG.csv", col_names = TRUE)
create_barplot(HOG)
create_boxplot(HOG, metadata)

cAMP <- read_csv("results/fungal_pop/cAMP.csv", col_names = TRUE)
create_barplot(cAMP)
create_boxplot(cAMP, metadata)

calcineurin <- read_csv("results/fungal_pop/calcineurin.csv", col_names = TRUE)
create_barplot(calcineurin)
create_boxplot(calcineurin, metadata)

#### ALL impact variants ####
# ALL <- read_csv("data/ALL_filtered.csv", col_names = TRUE)
ALL <- read_csv("data/effects_filtered_quality.csv", col_names = TRUE)
ALL_filtered <- ALL
metadata <- metadata %>%
    select(strain, lineage)

# # Filter out variants that are present in all strains of a lineage
# strains_per_lineage <- metadata %>%
#     group_by(lineage) %>%
#     summarise(total_strains = n_distinct(strain))

# strains_per_var <- ALL %>%
#     group_by(lineage, var_id, impact) %>%
#     summarise(num_strains = n_distinct(strain))

# ALL_joined <- ALL %>%
#     select(lineage, strain, var_id, gene_id, impact)%>%
#     left_join(strains_per_var, by = c("var_id", "lineage", "impact"))%>%
#     left_join(strains_per_lineage, by = "lineage")

# ALL_filtered <- ALL_joined %>%
#     filter(num_strains != total_strains)

# # Number of filtered out variants
# length(unique(ALL_joined$var_id)) - length(unique(ALL_filtered$var_id))

# Number of genes with variants 
ALL_filtered$impact <- factor(ALL_filtered$impact, levels = c("HIGH","MODERATE", "LOW"))

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

i_colors <- brewer.pal(3, "Set1")

plot <- ggplot(genes_per_strain, aes(x = lineage, y = num_genes)) +
    geom_quasirandom(aes(color = impact)) +
    geom_boxplot(alpha = 0) +
    facet_wrap(~impact, scales = "free", ncol = 1) +
    theme(panel.grid = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2),
          legend.position = "none") +
    labs(title = paste("Number genes with variants per strain by impact", sep = ""),
         x = "Lineage",
         y = "Number of genes") +
    scale_color_manual(values = i_colors) +
    scale_y_continuous(labels = scales::comma)
plot
ggsave("boxplot_genes_vars.png", plot, width = 9, height = 9)

# =============================================================================

#### finding essential genes ####
genes <- read_csv("genes.csv", col_names = TRUE)
genes$gene_id <- as.factor(genes$gene_id)
billmyre<- read_csv("results/fungal_pop/media-1.csv", col_names = TRUE)
billmyre <- billmyre %>%
    select(gene_id = Gene, essentiality = `Essentiality Classification`)
billmyre_essential <- billmyre %>%
    filter(essentiality == "ESS")

ALL <- read_csv("results/fungal_pop/ALL_filtered.csv", col_names = TRUE)

HIGH <- ALL_joined %>%
    filter(impact %in% c("MODERATE","HIGH"))
    
HIGH_num_strains <- HIGH %>%
    group_by(lineage,gene_id,var_id)%>%
    summarise(num_strains = n_distinct(strain))%>%
    ungroup()

HIGH_summary <- HIGH_num_strains%>%
    group_by(gene_id, lineage)%>%
    summarise(num_vars = n_distinct(var_id), num_strains = mean(num_strains))

genes_not_in_HIGH <- genes %>%
    anti_join(HIGH_summary, by = c("lineage","gene_id"))

VNI_essential <- genes_not_in_HIGH %>%
    filter(lineage == "VNI")
length(unique(VNI_essential$gene_id))
VNII_essential <- genes_not_in_HIGH %>%
    filter(lineage == "VNII")
length(unique(VNII_essential$gene_id))
VNBI_essential <- genes_not_in_HIGH %>%
    filter(lineage == "VNBI")
length(unique(VNBI_essential$gene_id))
VNBII_essential <- genes_not_in_HIGH %>%
    filter(lineage == "VNBII")
length(unique(VNBII_essential$gene_id))

billmyre_and_VNI<- inner_join(VNI_essential, billmyre_essential, by = "gene_id")
billmyre_and_VNII<- inner_join(VNII_essential, billmyre_essential, by = "gene_id")
billmyre_and_VNBI<- inner_join(VNBI_essential, billmyre_essential, by = "gene_id")
billmyre_and_VNBII<- inner_join(VNBII_essential, billmyre_essential, by = "gene_id")

#### Attempt to do Heatmap ####
samples <- read_csv("results/fungal_pop/samples.csv", col_names = TRUE)
metadata <- read_csv("/FastData/czirion/Crypto_Desjardins/config/metadata.csv", col_names = TRUE)
PWY <- read_csv("results/fungal_pop/HOG.csv", col_names = TRUE)
PWY_row_metadata <- metadata %>%
    select(sample, lineage) %>%
    distinct()%>%
    column_to_rownames("sample")
# PWY <- PWY %>%
#     filter(impact == "HIGH")

PWY_presence <- PWY %>%
    select(sample, var_id) %>%
    mutate(var_id = as.character(var_id)) %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = var_id, values_from = value, values_fill = 0)
# PWY_presence <- full_join(samples, PWY_presence, by = "sample") 
# PWY_presence[is.na(PWY_presence)] <- 0
PWY_presence <- PWY_presence %>%
    column_to_rownames("sample")
PWY_presence_matrix <- as.matrix(PWY_presence)

PWY_metadata <- PWY %>%
    select(var_id, gene_name, lineage) %>%
    distinct()%>%
    column_to_rownames("var_id")
PWY_row_metadata <- PWY %>%
    select(sample, lineage) %>%
    distinct()%>%
    column_to_rownames("sample")

col_ha <- columnAnnotation(Gene = PWY_metadata$gene_name)

Heatmap(PWY_presence_matrix,
    bottom_annotation = col_ha,
    column_split = PWY_metadata$gene_name,
    row_split = PWY_row_metadata$lineage)