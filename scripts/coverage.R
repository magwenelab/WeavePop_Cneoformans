log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)

print("Reading BED file")
all_chroms<- read.delim(snakemake@input[[1]], header = FALSE, col.names = c("Chromosome", "Start", "End", "Depth"), stringsAsFactors = TRUE)
#setwd("./genomes-annotations/SRS404449/")
#all_chroms<- read.delim("coverage_good.regions.bed.gz", header = FALSE, col.names = c("Chromosome", "Start", "End", "Depth"), stringsAsFactors = TRUE)

all_chroms <- all_chroms %>%
  group_by(Chromosome)%>%
  mutate(RoundDepth = round(Depth))%>%
  filter(RoundDepth != 0)

chrom_colors <- c(brewer.pal(nlevels(all_chroms$Chromosome)/2, "Dark2"),brewer.pal(nlevels(all_chroms$Chromosome)/2, "Set2"))

print("Plotting chromosome depth")
plot <- ggplot(data = all_chroms)+
  geom_col(aes(x= Start, y = RoundDepth, color = Chromosome, fill = Chromosome))+
  facet_wrap(~Chromosome,ncol = 2, scales = "free_x")+
  scale_color_manual(values = chrom_colors)+
  scale_fill_manual(values = chrom_colors)+
  scale_y_log10(name = "Coverage (X)", labels = comma)+
  scale_x_continuous(name = "Position (bp) ", labels = comma)+
  theme_light()+
  theme(legend.position="none")

#ggsave("coverage.svg", plot = plot, dpi = 200, units = "cm", height = 22, width = 22)

print("Saving plot")
ggsave(snakemake@output[[2]], plot = plot, dpi = 200, units = "cm", height = 22, width = 22)

print("Getting mean and median")
stats <- all_chroms %>%
  summarise(Mean = mean(Depth),
            Median = median(Depth))

stats <- stats %>%
pivot_longer(c(Mean, Median), names_to = "Stat", values_to = "Value")

plot <- ggplot(data = stats)+
  geom_point(aes(x = Chromosome, y = Value, color = Chromosome, shape = Stat))+
  scale_color_manual(values = chrom_colors)

ggsave(snakemake@output[[3]], plot = plot, dpi = 200, units = "cm", height = 10, width = 10)

write_csv(stats, snakemake@output[[1]])

print("Done!")