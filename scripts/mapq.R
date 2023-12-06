log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)
library(ggnewscale)

print("Reading BED file")


sample <- snakemake@input[[1]]
Split <- str_split(sample, "/")
sample <- Split[[1]][length(Split[[1]])-1]

raw<- read.delim(snakemake@input[[1]], header = FALSE, col.names = c("Accession", "Start", "End", "MAPQ"), stringsAsFactors = TRUE)
chrom_names <- read.csv(snakemake@input[[2]], header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
raw <- left_join(raw, chrom_names, by = "Accession")

lineage <- levels(as.factor(raw$Lineage))

raw_color = "aquamarine4"

print("Plotting chromosome depth")
plot <- ggplot()+
  geom_point(data = raw, aes(x= Start, y = MAPQ), size = 1, color = raw_color)+
  facet_wrap(~Chromosome,ncol = 2, scales = "free_x")+
  scale_x_continuous(name = "Position (bp) ", labels = comma)+
  theme_light()+
  theme(legend.position="right")+
  labs(title = paste(lineage, sample,  sep = " "))

ggsave(snakemake@output[[1]], plot = plot, dpi = 50, units = "cm", height = 22, width = 22)

