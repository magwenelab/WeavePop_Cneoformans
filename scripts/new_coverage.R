suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)
library(ggnewscale)

print("Reading BED file")

setwd("./genomes-annotations/SRS885880/")
raw<- read.delim("coverage.regions.bed.gz", header = FALSE, col.names = c("Accession", "Start", "End", "Depth"), stringsAsFactors = TRUE)
good<- read.delim("coverage_good.regions.bed.gz", header = FALSE, col.names = c("Accession", "Start", "End", "Depth"), stringsAsFactors = TRUE)
chrom_names <- read.csv("../../chromosome_names.csv", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
raw <- left_join(raw, chrom_names, by = "Accession")
good <- left_join(good, chrom_names, by = "Accession")

loci <- read.delim("../../results/loci_interest.tsv", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
loci <- loci %>% rename(Accession = seq_id)
loci <- left_join(loci, chrom_names, by = "Accession")

lineage <- levels(as.factor(raw$Lineage))
loci <- loci %>% filter(Lineage %in% lineage)

filtered_raw <- raw %>%
  group_by(Chromosome)%>%
  mutate(RoundDepth = round(Depth))%>%
  filter(RoundDepth != 0)

filtered_good <- good%>%
  group_by(Chromosome)%>%
  mutate(RoundDepth = round(Depth))%>%
  filter(RoundDepth != 0)

raw_color = "lightskyblue1"
good_color = "lightskyblue3"
color_quality = c("Good quality alignments" = good_color, "All alignments" = raw_color)

print("Plotting chromosome depth")
plot <- ggplot()+
  geom_col(data = filtered_raw, aes(x= Start, y = RoundDepth, fill= "All alignments" ))+ 
  geom_col(data = filtered_good, aes(x= Start, y = RoundDepth, fill = "Good quality alignments"))+
  scale_fill_manual(name= "Alignment quality", values= color_quality)+ 
  geom_point(data = loci, aes(x= start, y = 1, color = Loci), size = 1, shape = 15)+
  facet_wrap(~Chromosome,ncol = 2, scales = "free_x")+
  scale_y_log10(name = "Coverage (X)", labels = comma)+
  scale_x_continuous(name = "Position (bp) ", labels = comma)+
  theme_light()+
  theme(legend.position="right")

ggsave("../../coverage.svg", plot = plot, dpi = 50, units = "cm", height = 22, width = 22)

print("Getting mean and median")
stats_raw <- raw %>%
  group_by(Chromosome)%>%
  summarise(Mean = mean(Depth),
            Median = median(Depth))%>%
  pivot_longer(c(Mean, Median), names_to = "Measurement", values_to = "Value")

global_stats_raw<- raw %>%
  summarise(Mean = mean(Depth),
            Median = median(Depth))

stats_good <- good %>%
  group_by(Chromosome)%>%
  summarise(Mean = mean(Depth),
            Median = median(Depth))%>%
  pivot_longer(c(Mean, Median), names_to = "Measurement", values_to = "Value")

global_stats_good<- good %>%
  summarise(Mean = mean(Depth),
            Median = median(Depth))

toplim <- max(stats_raw$Value) + max(stats_raw$Value)/10

plot <- ggplot()+
  ylim(0,toplim) +
  geom_hline(aes(yintercept = global_stats_raw$Median,linetype = "Global median", color= "All alignments"))+
  geom_hline(aes(yintercept = global_stats_raw$Mean, linetype = "Global mean", color= "All alignments"))+
  geom_hline(aes(yintercept = global_stats_good$Median,linetype = "Global median", color = "Good quality alignments"))+
  geom_hline(aes(yintercept = global_stats_good$Mean, linetype = "Global mean", color = "Good quality alignments"))+
  geom_point(data = stats_raw, aes(x = factor(Chromosome, levels = as.character(sort(unique(Chromosome)))), y = Value, shape = Measurement, color= "All alignments"))+ 
  geom_point(data = stats_good, aes(x = factor(Chromosome, levels = as.character(sort(unique(Chromosome)))), y = Value, shape = Measurement, color = "Good quality alignments"))+ 
  labs(y = "Coverage", x = "Chromosome")+
  theme_light()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_shape_manual(values = c(16,15), name = NULL)+
  scale_linetype_manual(values = c("solid","dotted"), name = NULL)+
  scale_color_manual(name= "Alignment quality", values= color_quality)

ggsave("../../coverage_stats.svg", plot = plot, dpi = 50, units = "cm", height = 15, width = 15)

write_csv(stats_raw, "../../coverage.txt")

print("Done!")