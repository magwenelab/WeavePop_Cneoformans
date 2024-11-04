suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggbeeswarm))


setwd("/FastData/czirion/Crypto_Desjardins/fungal_pop/ashton_desj")


metadata <- read_csv("metadata.csv", col_names = TRUE)
metadata <- metadata %>%
    select(strain, lineage)

strains_per_lineage <- metadata %>%
    group_by(lineage) %>%
    summarise(total_strains = n_distinct(strain))

HML <- read_csv("HML.csv", col_names = TRUE)
HML$impact <- factor(HML$impact, levels = c("HIGH", "MODERATE", "LOW"))

strains_per_var <- HML %>%
    group_by(lineage, var_id, impact) %>%
    summarise(num_strains = n_distinct(strain))

HML_joined <- HML %>%
    select(lineage, strain, var_id, gene_id, impact)%>%
    left_join(strains_per_var, by = c("var_id", "lineage", "impact"))%>%
    left_join(strains_per_lineage, by = "lineage")

HML_filtered <- HML_joined %>%
    filter(num_strains != total_strains)


number_genes_vars <- HML_filtered %>%
    group_by(strain, impact, gene_id) %>%
    summarise(num_vars = n_distinct(var_id)) %>%
    ungroup() %>%
    group_by(strain, impact) %>%
    summarise(num_genes = n_distinct(gene_id), mean_vars = mean(num_vars))

expanded <- expand.grid(
    strain = unique(metadata$strain),
    impact = unique(HML_filtered$impact))

genes_per_strain <- left_join(expanded, number_genes_vars, by = c("strain", "impact"))
genes_per_strain$num_genes[is.na(genes_per_strain$num_genes)] <- 0

genes_per_strain <- left_join(genes_per_strain, metadata, by = "strain")

i_colors <- brewer.pal(3, "Set1")

plot <- ggplot(genes_per_strain, aes(x = lineage, y = num_genes, alpha=mean_vars)) +
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
ggsave("HML_boxplot_genes.png", plot, width = 9, height = 9)
