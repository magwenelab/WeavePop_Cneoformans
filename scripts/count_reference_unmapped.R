log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(pheatmap)
library(RColorBrewer)

unmapped<-read.csv(snakemake@input[[1]], header = FALSE, col.names = c("gene", "lineage"))
unmapped_count <- unmapped %>% 
    group_by(lineage) %>% 
    count()
write.csv(unmapped_count, snakemake@output[[1]], row.names = FALSE, quote = FALSE)

#unmapped<-read.csv("reference_genomes_65/references_unmapped_features.csv", header = FALSE, col.names = c("gene", "lineage"))
#unmapped_count <- unmapped %>% 
#    group_by(lineage) %>% 
#    count()
#write.csv(unmapped_count, "reference_genomes_65/references_unmapped_count.csv", row.names = FALSE)

genes<-read_delim(snakemake@input[[2]], col_names = TRUE, na = "N/A", show_col_types = FALSE )
genes <- genes %>% as.data.frame() %>% select(gene = ID, Name, description, primary_tag, seq_id)

rownames(genes)<- genes$gene
genes<- genes %>% select(primary_tag, seq_id)

unmapped$unmapped <- 0
unmapped_wide <- as.data.frame(unmapped %>% pivot_wider(names_from = lineage, values_from = unmapped, values_fill = 1))

#unmapped_wide <- left_join(unmapped_wide, genes, by = "gene")

rownames(unmapped_wide)<- unmapped_wide$gene
unmapped_wide<- unmapped_wide %>% select(-gene)
unmapped_matrix<-as.matrix(unmapped_wide)

lg.brks<-c(0,1) 
lb.brks<-c("Unmapped","Mapped") 
plot <- pheatmap(unmapped_matrix,
            color =  c("lightgray", "hotpink4"),
            annotation_row = genes,
            legend_breaks = lg.brks, 
            legend_labels = lb.brks,
            clustering_method = "median",
            cellwidth = 10,
            cellheight = 2,
            fontsize_row= 2,
            treeheight_row = 0,
            treeheight_col = 0,
            filename = snakemake@output[[2]])

#ggsave("unmapped_references.svg", plot, dpi = 200, units = "cm", height = 23, width = 16)
#ggsave(snakemake@output[[2]], plot, dpi = 200, units = "cm", height = 23, width = 16)

# Missing to add gene name labels