log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)

print("Reading BED file")

setwd("./genomes-annotations/SRS885169/")
raw<- read.delim("coverage.regions.bed.gz", header = FALSE, col.names = c("Accession", "Start", "End", "Depth"), stringsAsFactors = TRUE)
good<- read.delim("coverage_good.regions.bed.gz", header = FALSE, col.names = c("Accession", "Start", "End", "Depth"), stringsAsFactors = TRUE)

chrom_names <- read.csv("../../chromosome_names.csv", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
loci <- read.delim("../../results/VNBI_loci_interest.tsv", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
loci <- loci %>% rename(Accession = seq_id)
loci <- left_join(loci, chrom_names, by = "Accession")

raw <- left_join(raw, chrom_names, by = "Accession")
good <- left_join(good, chrom_names, by = "Accession")

filtered_raw <- raw %>%
  group_by(Chromosome)%>%
  mutate(RoundDepth = round(Depth))%>%
  filter(RoundDepth != 0)

filtered_good <- good%>%
  group_by(Chromosome)%>%
  mutate(RoundDepth = round(Depth))%>%
  filter(RoundDepth != 0)

print("Plotting chromosome depth")
plot <- ggplot()+
  geom_col(data = filtered_raw, aes(x= Start, y = RoundDepth), color = "lightskyblue1")+ 
  geom_col(data = filtered_good, aes(x= Start, y = RoundDepth), color = "lightskyblue3")+ 
  geom_point(data = loci, aes(x= start, y = 1), size = 1, shape = 15)+
  facet_wrap(~Chromosome,ncol = 2, scales = "free_x")+
  scale_y_log10(name = "Coverage (X)", labels = comma)+
  scale_x_continuous(name = "Position (bp) ", labels = comma)+
  theme_light()+
  theme(legend.position="right")

ggsave("../../coverage.svg", plot = plot, dpi = 200, units = "cm", height = 22, width = 22)

print("Getting mean and median")
stats <- raw %>%
  group_by(Chromosome)%>%
  summarise(Mean = mean(Depth),
            Median = median(Depth))

stats <- stats %>%
pivot_longer(c(Mean, Median), names_to = "Measurement", values_to = "Value")
global_stats<- raw %>%
  summarise(Mean = mean(Depth),
            Median = median(Depth))

toplim <- max(stats$Value) + max(stats$Value)/10

plot <- ggplot(data = stats, aes(x = Chromosome, y = Value, shape = Measurement))+
  ylim(0,toplim) +
  geom_hline(aes(yintercept = global_stats$Median,linetype = "Global median"))+
  geom_hline(aes(yintercept = global_stats$Mean, linetype = "Global mean"))+
  geom_point()+ 
  labs(y = "Coverage")+
  theme_light()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_shape_manual(values = c(16,15), name = NULL)+
  scale_linetype_manual(values = c("solid","dotted"), name = NULL)

ggsave("../../coverage_stats.svg", plot = plot, dpi = 200, units = "cm", height = 10, width = 10)

write_csv(stats, "coverage.txt")

print("Done!")