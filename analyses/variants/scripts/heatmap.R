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

get_matrix_filter_genes <- function(vars_df, metadata){

    matrix <- vars_df%>%
        select(gene_id, strain, presence_vars)%>%
        pivot_wider(names_from = strain, values_from = presence_vars) %>%
        column_to_rownames("gene_id")%>%
        as.matrix()

    matrix <- matrix[rowSums(matrix) > 0,]
    matrix <- matrix[!is.na(rownames(matrix)), ]
    
    vars_df <- vars_df%>%
        filter(gene_id %in% rownames(matrix))
    
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

    colors_df <- metadata %>%
        select(lineage, color)%>%
        distinct()%>%
        arrange(lineage)
    lineage_colors <- colors_df$color
    names(lineage_colors) <- colors_df$lineage

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
    
    colors_ess_df <- list_data$df %>%
        select(essentiality, color)%>%
        distinct()%>%
        arrange(essentiality)
    ess_colors <- colors_ess_df$color
    names(ess_colors) <- colors_ess_df$essentiality

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

plot_boxplot <- function(data, filename){

    genes_per_strain <- data %>%
        group_by(lineage, strain) %>%
        summarise(num_genes = sum(presence_vars))

    boxplot <- ggplot(genes_per_strain, aes(x = lineage, y = num_genes)) +
        geom_boxplot() +
        geom_quasirandom(aes(color = lineage)) +
        facet_wrap(~lineage, nrow  = 1, scales = "free_x")+
        theme(panel.grid = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2),
            legend.position = "none")+
        labs(title = "Number of genes with variants per strain",
                x = "Lineage",
                y = "Number of genes")
    ggsave(filename, boxplot, height = 10, width = 10, units = "in", dpi = 300)

}

plot_tree <- function(data, metadata, tree_file, filename){
    metadata <- metadata %>%
        select(strain, lineage)%>%
        add_row(strain = "H99", lineage = "VNI")
        
    tree_des <- read.tree(tree_file)
    outgroup_clade <- metadata %>%
        filter(lineage == "VNII")%>%
        select(strain)
    outgroup_clade <- outgroup_clade$strain
    tree_des <- root(tree_des, outgroup = outgroup_clade)

    genes_per_strain <- data %>%
        select(strain, presence_vars)%>%
        group_by(strain) %>%
        summarise(num_genes = sum(presence_vars)) 

    plot_tree <- ggtree(tree_des) %<+% metadata + 
                    geom_tiplab(aes(label = label), 
                    size = 2, align = TRUE)

    var_tree <- plot_tree + 
                geom_facet(panel= "Number of genes with variants",
                            data = genes_per_strain,
                            geom = geom_col,
                            mapping = aes(x = num_genes, color = lineage, fill = lineage),
                            orientation = "y")+
                theme_tree2(legend.position=c(.05, .85))+
                theme(panel.border = element_blank())+
                guides(fill = guide_legend(title = "Lineage"), color = guide_legend(title = "Lineage"))
    ggsave(filename, var_tree, height = 20, width = 10, units = "in", dpi = 300)

}

#### Input ####
setwd("/FastData/czirion/Crypto_Diversity_Pipeline/")

metadata <- read_csv("Crypto_Desjardins/config/metadata.csv", col_names = TRUE) 
metadata <- metadata %>%
    mutate(color = brewer.pal(n = 8, name = "Set2")[match(lineage, unique(lineage))])%>%
    as.data.frame()

billmyre<- read_csv("analyses/variants/data/media-1.csv", col_names = TRUE)
billmyre <- billmyre %>%
    select(gene_id = Gene, essentiality = `Essentiality Classification`) 
billmyre <- billmyre %>%
        mutate(essentiality = as.factor(essentiality))%>%
        mutate(essentiality = recode(essentiality, "NESS" = "Non Essential",
                                                    "ESS" = "Essential",
                                                    "UNK" = "Unknown"))%>%
        mutate(color = rev(brewer.pal(n = 12, name = "Paired"))[match(essentiality, levels(essentiality))])%>%
        as.data.frame()
                                    
                                            
sets <- list.files(path = "analyses/variants/data/variants_per_gene_per_strain", pattern = "*.csv", full.names = TRUE)

for (set in sets){
    print(paste("Analizing set:", set))
    set_name <- gsub("analyses/variants/data/variants_per_gene_per_strain", "", set)
    set_name <- gsub(".csv", "", set_name)
    heatmap_filename <- paste("analyses/variants/results_heatmaps/", set_name, ".png", sep = "")
    boxplot_filename <- paste("analyses/variants/results_boxplots/", set_name, ".png", sep = "")
    tree_barplot_filename <- paste("analyses/variants/results_tree_barplots/", set_name, ".png", sep = "")

    if (!file.exists(heatmap_filename)) {
        print("Plotting heatmap")
        presence <- read_csv(set, col_names = TRUE, show_col_types = FALSE)
        presence <- left_join(presence, billmyre, by = "gene_id")
        data <-  get_matrix_filter_genes(presence, metadata)
        ann <- get_heatmap_annotation(data)
        plot_heatmap(data, ann, heatmap_filename)
    } else {
        print("Heatmap already plotted")
    }

    if (!file.exists(boxplot_filename)) {
        print("Plotting boxplot")
        presence <- read_csv(set, col_names = TRUE, show_col_types = FALSE)
        plot_boxplot(presence, boxplot_filename)
    } else {
        print("Boxplot already plotted")
    }

    if (!file.exists(tree_barplot_filename)) {
        print("Plotting tree barplot")
        presence <- read_csv(set, col_names = TRUE, show_col_types = FALSE)
        plot_tree(presence, metadata, "analyses/variants/data/CryptoDiversity_Desjardins_Tree.tre", tree_barplot_filename)
    } else {
        print("Tree barplot already plotted")
    }

}

