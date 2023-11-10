log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)

print("Reading BED file")
all_chroms<- read.delim("Downloads/SRS404442.regions.bed", header = FALSE, col.names = c("Chromosome", "Start", "End", "Depth"), stringsAsFactors = TRUE)

all_chroms<- read.delim(snakemake@input[[1]], header = FALSE, col.names = c("Chromosome", "Start", "End", "Depth"), stringsAsFactors = TRUE)

all_chroms <- all_chroms %>%
  group_by(Chromosome)

print("Making file of average depth")
summary <- all_chroms %>%
  summarise(Avg_Depth = mean(Depth))

write_csv(summary, snakemake@output[[1]])

chrom_colors <- c( brewer.pal(nlevels(all_chroms$Chromosome)/2, "Dark2"),brewer.pal(nlevels(all_chroms$Chromosome)/2, "Set2"))

print("Plotting chromosome depth")
plot <- ggplot(data = all_chroms)+
  geom_col(aes(x= Start, y = Depth, color = Chromosome, fill = Chromosome))+
  facet_wrap(~Chromosome,ncol = 2, scales = "free_x")+
  scale_color_manual(values = chrom_colors)+
  scale_fill_manual(values = chrom_colors)+
  scale_y_log10(name = "Depth (X)")+
  scale_x_continuous(name = "Position (bp) ", labels = comma)+
  theme_light()+
  theme(legend.position="none")

print("Saving plot")
ggsave(snakemake@output[[2]], plot = plot, dpi = 200, units = "cm", height = 22, width = 22)
print("Done!")