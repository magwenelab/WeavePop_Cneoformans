log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)

print("Reading TSV file")
mapq<- read.csv(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE)
#setwd("./genomes-annotations/SRS881221/")
#mapq<- read.csv("mapq.csv", header = TRUE, stringsAsFactors = TRUE)

mapq <-mapq %>%
  rename(Accession = Chromosome)
#chrom_names <- read.csv("../../chromosome_names.csv", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))

chrom_names <- read.csv("chromosome_names.csv", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))

mapq <- left_join(mapq, chrom_names, by = "Accession")


print("Plotting MAPQ distribution")

#chrom_colors <- c(brewer.pal(nlevels(mapq$Chromosome)/2, "Dark2"),brewer.pal(nlevels(mapq$Chromosome)/2, "Set2"))

plot <- ggplot(mapq, aes(x=MAPQ, y=Count))+
  geom_col(aes(fill = Chromosome))+
  #scale_fill_manual(values = chrom_colors)+
  facet_wrap(~Chromosome,ncol = 2)+
  scale_y_log10(name = "Number of Reads", labels = comma, breaks = 100 * 10^seq(0,4, by = 2))+
  scale_x_continuous(name = "Mapping Quality", n.breaks = 8)+
  theme_light()+
  theme(legend.position="none")

#ggsave("mapq_distribution.svg", plot = plot, dpi = 200, units = "cm", height = 22, width = 22)

print("Saving plot")
ggsave(snakemake@output[[1]], plot = plot, dpi = 200, units = "cm", height = 22, width = 22)
print("Done!")

# Binned version

#mapq_bin <- mapq %>% 
#  mutate(bin_mapq = ntile(MAPQ,n=6))%>%
#  group_by(Chromosome, bin_mapq)%>%
#  reframe(binned_count = sum(Count))

#plot <- ggplot(mapq_bin, aes(x=bin_mapq, y=binned_count))+
#  geom_col(aes(fill = Chromosome))+
#  #scale_fill_manual(values = chrom_colors)+
#  facet_wrap(~Chromosome,ncol = 2)+
#  scale_y_log10(name = "Number of Reads", labels = comma, breaks = 100 * 10^seq(0,4, by = 2))+
#  scale_x_continuous(name = "Mapping Quality", n.breaks = 8)+
#  theme_light()+
#  theme(legend.position="none")

##ggsave("mapq_distribution_binned.svg", plot = plot, dpi = 200, units = "cm", height = 22, width = 22)