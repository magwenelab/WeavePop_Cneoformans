log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)

print("Reading TSV file")
mapq<- read.delim(snakemake@input[[1]], header = FALSE, col.names = c("Chromosome", "MAPQ", "Count"), stringsAsFactors = TRUE)

print("Plotting MAPQ distribution")

chrom_colors <- c(brewer.pal(nlevels(mapq$Chromosome)/2, "Dark2"),brewer.pal(nlevels(mapq$Chromosome)/2, "Set2"))

plot <- ggplot(mapq, aes(x=MAPQ, y=Count))+
  geom_col(aes( fill = Chromosome))+
  scale_fill_manual(values = chrom_colors)+
  facet_wrap(~Chromosome,ncol = 2)+
  scale_y_continuous(name = "Number of Reads", labels = comma)+
  scale_x_continuous(name = "Mapping Quality")+
  theme_light()+
  theme(legend.position="none")

print("Saving plot")
ggsave(snakemake@output[[2]], plot = plot, dpi = 200, units = "cm", height = 22, width = 22)
print("Done!")