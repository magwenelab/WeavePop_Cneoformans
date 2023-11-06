log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)

print("Reading TSV file")
cov<- read.delim("Downloads/SRS404442/cov.tsv", header = FALSE, col.names = c("Chromosome", "Range", "Coverage", "Count"), stringsAsFactors = TRUE)

print("Plotting Coverage distribution")

chrom_colors <- c(brewer.pal(nlevels(cov$Chromosome)/2, "Dark2"),brewer.pal(nlevels(cov$Chromosome)/2, "Set2"))

plot <- ggplot(cov, aes(x=Coverage, y=Count))+
  geom_line(aes(color = Chromosome))+
  scale_color_manual(values = chrom_colors)+
  facet_wrap(~Chromosome,ncol = 2)+
  scale_y_continuous(name = "Number of Reads", labels = comma)+
  theme_light()+
  theme(legend.position="none")

print("Saving plot")
ggsave(snakemake@output[[1]], plot = plot, dpi = 200, units = "cm", height = 22, width = 22)
print("Done!")