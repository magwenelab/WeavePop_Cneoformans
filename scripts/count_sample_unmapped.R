log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(pheatmap)

unmapped<-read.csv(snakemake@input[[1]], header = FALSE, col.names = c("gene", "sample"))
unmapped_count <- unmapped %>% 
    group_by(sample) %>% 
    count()
write.csv(unmapped_count, snakemake@output[[1]], row.names = FALSE, quote = FALSE)

#unmapped<-read.csv("samples_unmapped_features.csv", header = FALSE, col.names = c("gene", "sample"))
#unmapped_count <- unmapped %>% 
#   group_by(sample) %>% 
#    count()
#write.csv(unmapped_count, "reference_genomes_65/references_unmapped_count.csv", row.names = FALSE)

genes<-read_delim(snakemake@input[[2]], col_names = TRUE, na = "N/A", show_col_types = FALSE )
genes <- genes %>% as.data.frame() %>% select(gene = ID, Name, description, primary_tag, seq_id)
rownames(genes)<- genes$gene
genes<- genes %>% select(primary_tag, seq_id)

metadata<-read_csv(snakemake@input[[3]], col_names = TRUE, show_col_types = FALSE)
metadata <- metadata %>% as.data.frame() %>% select(sample, lineage)
rownames(metadata)<- metadata$sample
metadata<- metadata %>% select(lineage)

unmapped$unmapped <- 0
unmapped_wide <- as.data.frame(unmapped %>% pivot_wider(names_from = sample, values_from = unmapped, values_fill = 1))
rownames(unmapped_wide)<- unmapped_wide$gene
unmapped_wide<- unmapped_wide %>% select(-gene)
unmapped_matrix<-as.matrix(unmapped_wide)
lg.brks<-c(0,1) 
lb.brks<-c("Unmapped","Mapped") 
plot <- pheatmap(unmapped_matrix,
            color =  c("lightgray", "hotpink4"),
            annotation_row = genes,
            annotation_col = metadata,
            legend_breaks = lg.brks, 
            legend_labels = lb.brks,
            clustering_method = "median",
            cellwidth = 10,
            cellheight = 2,
            fontsize_row= 2,
            treeheight_row = 0,
            treeheight_col = 0,
            filename = snakemake@output[[2]])

#ggsave("unmapped_references.svg", plot, dpi = 200, units = "cm", height = 22, width = 10)
#ggsave(snakemake@output[[2]], plot, dpi = 200, units = "cm", height = 22, width = 10)
