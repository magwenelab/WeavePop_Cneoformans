suppressPackageStartupMessages(library(tidyverse))

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