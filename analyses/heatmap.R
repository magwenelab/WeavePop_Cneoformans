suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))

#### Get matrix of number of variants in each gene in each strain ####

get_matrix <- function(vars_df, metadata){

    matrix <- vars_df%>%
        select(gene_id, strain, presence_vars)%>%
        pivot_wider(names_from = strain, values_from = presence_vars) %>%
        column_to_rownames("gene_id")%>%
        as.matrix()
    
    metadata <- metadata %>%
        filter(strain %in% colnames(matrix))%>%
        mutate(strain = factor(strain, levels = colnames(matrix))) %>%
        arrange(strain)

    list <- list(matrix = matrix, df = vars_df, metadata = metadata)
    return(list)
}

#### Get annotation for heatmap ####
get_heatmap_annotation <- function(list_data){

    ann <- data.frame(list_data$metadata$lineage)
    colnames(ann) <- 'Lineage'
    num_lineages <- length(unique(list_data$metadata$lineage))
    lineage_colors <- brewer.pal(n = num_lineages, name = "Set2")#[1:num_lineages]
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
    
    num_ess <- length(unique(ess$essentiality))
    ess_colors <- rev(brewer.pal(n = 12, name = "Paired"))[1:num_ess]
    names(ess_colors) <- unique(ess$essentiality)
    print(ess_colors)

    ess <- data.frame(ess$essentiality)
    colnames(ess)<- 'Essentiality'

    row_ha <- rowAnnotation(df = ess,
                col = list(Essentiality = ess_colors))

    col_fun <- colorRamp2(c(0,1), c("gray95", "darkblue"))

    list <- list(col_ha = col_ha, row_ha = row_ha, col_fun = col_fun)
    return(list)
}

#### Plot heatmap ####
plot_heatmap <- function(list_data, list_ann, filename){
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
                heatmap_legend_param = list(title = "Presence of Variants", na_col = "black"),
                row_title = "Genes")

    png(file = filename, width = 16, height = 9, units = "in", res = 1200)
    draw(plot,
        annotation_legend_side = "bottom")
    dev.off()
}


#### Input ####
setwd("/FastData/czirion/Crypto_Desjardins/fungal_pop")

metadata <- read_csv("../config/metadata.csv", col_names = TRUE)
billmyre<- read_csv("data/media-1.csv", col_names = TRUE)
billmyre <- billmyre %>%
    select(gene_id = Gene, essentiality = `Essentiality Classification`) 
billmyre <- billmyre %>%
        mutate(essentiality = as.factor(essentiality))%>%
        mutate(essentiality = recode(essentiality, "ESS" = "Essential",
                                                    "NESS" = "Non Essential",
                                                    "UNK" = "Unknown"))
sets <- list.files(path = "data/variants_per_gene", pattern = "*.csv", full.names = TRUE)

for (set in sets){
    print(paste("Analizing set:", set))
    set_name <- gsub("data/variants_per_gene/", "", set)
    set_name <- gsub(".csv", "", set_name)
    plot_filename <- paste("plots/", set_name, ".png", sep = "")
    print(plot_filename)
    
    if (!file.exists(plot_filename)) {
        print("Plotting heatmap")
        presence <- read_csv(set, col_names = TRUE, show_col_types = FALSE)
        presence <- left_join(presence, billmyre, by = "gene_id")
        data <-  get_matrix(presence, metadata)
        ann <- get_heatmap_annotation(data)
        plot_heatmap(data, ann, plot_filename)
    } else {
        print("Heatmap already plotted")
    }
}
