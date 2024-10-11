suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))

setwd("/FastData/czirion/Crypto_Desjardins/fungal_pop")

#### Input ####
metadata <- read_csv("../config/metadata.csv", col_names = TRUE, show_col_types = FALSE)
annotation <- read_csv("data/genes.csv", col_names = TRUE, show_col_types = FALSE)
effects <- read_csv("data/ALL_filtered.csv", col_names = TRUE, show_col_types = FALSE)
billmyre<- read_csv("data/media-1.csv", col_names = TRUE)

#### Organize essential genes ####
billmyre <- billmyre %>%
    select(gene_id = Gene, essentiality = `Essentiality Classification`) 
billmyre <- billmyre %>%
        mutate(essentiality = as.factor(essentiality))%>%
        mutate(essentiality = recode(essentiality, "ESS" = "Essential",
                                                    "NESS" = "Non Essential",
                                                    "UNK" = "Unknown"))

#### Filter out variants that are present in all strains of a lineage ####
metadata <- metadata %>%
    select(strain, lineage)

strains_per_lineage <- metadata %>%
    group_by(lineage) %>%
    summarise(total_strains = n_distinct(strain))

strains_per_var <- effects %>%
    group_by(lineage,var_id) %>%
    summarise(num_strains = n_distinct(strain))

variants_filtered <- left_join(strains_per_var, strains_per_lineage, by = "lineage") %>%
        filter(num_strains != total_strains)

effects_filtered <- effects %>%
    filter(var_id %in% variants_filtered$var_id)

#### Get number of genes with variants ####
effects_filtered$impact <- factor(effects_filtered$impact,  levels = c("HIGH","MODERATE", "LOW"))

vars_per_gene <- effects_filtered %>%
    group_by(lineage, strain, impact, gene_id) %>%
    summarise(num_vars = n_distinct(var_id)) %>%
    ungroup()

#### Get list of genes present in all strains ####
genes_in_strains <- full_join(metadata, annotation, by = "lineage")

#### Get all combinations of genes, impacts and strains ####
genes_in_strains_impacts <- genes_in_strains %>%
    expand(nesting(lineage, gene_id, strain), impact = c("HIGH", "MODERATE", "LOW"))

#### Get number of variants in all combinations of genes, impacts and strains ####
#                   and number of strains with variants in gene

all_genes_vars <- left_join(genes_in_strains_impacts, vars_per_gene,
                    by = c("lineage", "strain", "gene_id", "impact"))%>%
    mutate(num_vars = replace_na(num_vars, 0))%>%
    mutate(any_variant = ifelse(num_vars == 0, 0, 1))%>%
    group_by(gene_id, impact, lineage)%>%
    mutate(strains_with_vars_in_gene = sum(any_variant))%>%
    ungroup()%>%
    left_join(strains_per_lineage, by = "lineage")%>%
    mutate(percent_strains_with_vars = strains_with_vars_in_gene / total_strains * 100 )%>%
    left_join(billmyre, by = "gene_id")

#### Get matrix of number of variants in each gene in each strain ####

get_matrix <- function(vars_df, 
                        impact_filter, lineage_filter, percent_strains_filter){
    
    df <- vars_df %>%
        filter(impact == impact_filter)%>%
        filter(lineage %in% lineage_filter)%>%
        filter(percent_strains_with_vars >= percent_strains_filter)%>%
        droplevels()

    matrix <- df %>%
        select(gene_id, strain, num_vars)%>%
        pivot_wider(names_from = strain, values_from = num_vars)%>%
        column_to_rownames("gene_id")%>%
        as.matrix()
    
    metadata <- metadata %>%
        filter(strain %in% colnames(matrix))%>%
        mutate(strain = factor(strain, levels = colnames(matrix))) %>%
        arrange(strain)

    list <- list(matrix = matrix, df = df, metadata = metadata)
    return(list)
}

get_heatmap_annotation <- function(list_data){

    ann <- data.frame(list_data$metadata$lineage)
    colnames(ann) <- 'Lineage'
    num_lineages <- length(unique(list_data$metadata$lineage))
    lineage_colors <- rev(brewer.pal(n = 8, name = "Set2"))[1:num_lineages]
    names(lineage_colors) <- unique(list_data$metadata$lineage)

    col_ha <- columnAnnotation(
                df = ann,
                col = list(Lineage = lineage_colors),
                annotation_legend_param = list(
                    Lineage = list(direction = "horizontal", nrow = 1)))

    ess <- list_data$df %>%
        select(gene_id, essentiality)%>%
        distinct()%>%
        mutate(gene_id = factor(gene_id, levels = rownames(list_data$matrix)))%>%
        arrange(gene_id)
    ess <- data.frame(ess$essentiality)
    colnames(ess)<- 'Essentiality'
    ess_colors <- brewer.pal(n = 3, name = "Dark2")
    names(ess_colors) <- unique(list_data$df$essentiality)

    row_ha <- rowAnnotation(df = ess,
                col = list(Essentiality = ess_colors))

    max_vars <- max(list_data$df$num_vars)
    palette <- brewer.pal(n = 5, name = "RdYlBu")
    col_fun <- colorRamp2(c(0,1,max_vars/4, max_vars/2, max_vars*3/4, max_vars ), c("gray95", palette))

    list <- list(col_ha = col_ha, row_ha = row_ha, col_fun = col_fun)
    return(list)
}

plot_heatmap <- function(list_data, list_ann){
    plot <- Heatmap(list_data$matrix,
                top_annotation = list_ann$col_ha,
                bottom_annotation = list_ann$col_ha,
                left_annotation = list_ann$row_ha,
                right_annotation = list_ann$row_ha,
                na_col = "black",
                col = list_ann$col_fun,
                show_row_names = FALSE,
                use_raster = TRUE,
                column_names_gp = gpar(fontsize = 3),
                heatmap_legend_param = list(title = "Number of Variants", na_col = "black"),
                column_title = "Number of variants per gene in each strain",
                row_title = "Genes")

    png(file = "heatmap.png", width = 16, height = 9, units = "in", res = 1200)
    draw(plot,
        annotation_legend_side = "bottom")
    dev.off()
}

HIGH_VNI <- get_matrix(all_genes_vars, "HIGH", c("VNII"), 10)
HIGH_VNI_ann <- get_heatmap_annotation(HIGH_VNI)
plot_heatmap(HIGH_VNI, HIGH_VNI_ann)



###################


#### Simple heatmap for testing ####


set.seed(123)  # For reproducibility

# Create a 10x10 matrix with random numbers
random_matrix <- matrix(sample(1:10, 10, replace = TRUE), nrow = 10, ncol = 10)

na_indices <- sample(length(random_matrix), 10)  # Randomly select 10 indices to be NA
zero_indices <- sample(length(random_matrix), 10)
random_matrix[na_indices] <- NA
random_matrix[zero_indices] <- 0

# Assign random letters as column names and row names
colnames(random_matrix) <- sample(LETTERS, 10)
rownames(random_matrix) <- sample(letters, 10)

print(random_matrix)

df <- as.data.frame(random_matrix)
df$sum_row <- rowSums(df, na.rm=TRUE)
df_1 <- df
df_1[df_1 != 0] <- 1
df_1$sum_row <- rowSums(df_1, na.rm=TRUE)

metad <- data.frame(strain = colnames(random_matrix),
    lineage = c("VNI", "VNI", "VNI", "VNB","VNB", "VNII", "VNII", "VNII","VNBB", "VNBB"))

metad <- metad %>%
    column_to_rownames("strain")

col_fun <- colorRamp2(c(0,1), c("white", "purple"))

lineage_colors <- brewer.pal(n = 4, name = "Dark2")
names(lineage_colors) <- unique(metad$lineage)

col_ha <- columnAnnotation(lineage = metad$lineage,
            col = list(lineage = lineage_colors),
            annotation_legend_param = list(
                lineage = list(direction = "horizontal", nrow = 1)))


ann <- data.frame(metad$lineage)
colnames(ann) <- 'Lineage'
col_ha <- columnAnnotation(df = ann)
lineage_colors <- brewer.pal(n = 4, name = "Dark2")
names(lineage_colors) <- unique(metad$lineage)

col_ha <- columnAnnotation(df = ann,
            col = list(Lineage = lineage_colors),
            annotation_legend_param = list(
                Lineage = list(direction = "horizontal", nrow = 1)))

plot <-Heatmap(random_matrix, 
    top_annotation = col_ha,
    col = col_fun,
    na_col = "black",
    heatmap_legend_param = list(title = "Number of Variants", na_col = "black"),
    column_title = "Number of high impact variants per gene in each strain")

draw(plot)

draw(plot,
    annotation_legend_side = "bottom")

