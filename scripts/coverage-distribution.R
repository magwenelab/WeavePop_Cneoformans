log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)

print("Reading TSV file")
cov<- read.csv(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE)
#setwd("./genomes-annotations/SRS404449/")
#cov<- read.csv("cov.csv", header = TRUE, stringsAsFactors = TRUE)

cov <-cov %>%
  rename(Accession = Chromosome)

#chrom_names <- read.csv("../../chromosome_names.csv", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
chrom_names <- read.csv("chromosome_names.csv", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
cov <- left_join(cov, chrom_names, by = "Accession")

print("Plotting Coverage distribution")

#chrom_colors <- c(brewer.pal(nlevels(cov$Chromosome)/2, "Dark2"),brewer.pal(nlevels(cov$Chromosome)/2, "Set2"))

plot <- ggplot(cov, aes(x=Coverage, y=Count))+
  geom_col(aes(color = Chromosome))+ 
  #scale_color_manual(values = chrom_colors)+
  facet_wrap(~Chromosome,ncol = 2)+
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Coverage (X) ", labels = comma, n.breaks = 10)+
  theme_light()+
  theme(legend.position="none")
#ggsave("cov_distribution.svg", plot = plot, dpi = 200, units = "cm", height = 22, width = 22)

print("Saving plot")
ggsave(snakemake@output[[1]], plot = plot, dpi = 200, units = "cm", height = 22, width = 22)
print("Done!")
