library(tidyverse)
library(RColorBrewer)
library(scales)
library(svglite)
all_chroms<- read.delim("Downloads/SRS404442.regions.bed", header = FALSE, col.names = c("Chromosome", "Start", "End", "Depth"), stringsAsFactors = TRUE)

all_chroms <- all_chroms %>%
  group_by(Chromosome)

summary <- all_chroms %>%
  summarise(Avg.Depth = mean(Depth))

write_csv()

chrom_colors <- c( brewer.pal(nlevels(all_chroms$Chromosome)/2, "Dark2"),brewer.pal(nlevels(all_chroms$Chromosome)/2, "Set2"))

plot <- ggplot(data = all_chroms)+
  geom_col(aes(x= Start, y = Depth, color = Chromosome, fill = Chromosome))+
  facet_wrap(~Chromosome,ncol = 2, scales = "free_x")+
  scale_color_manual(values = chrom_colors)+
  scale_fill_manual(values = chrom_colors)+
  scale_y_continuous(name = "Depth (X)", labels = comma)+
  scale_x_continuous(name = "Position (bp) ", labels = comma)+
  theme_minimal()+
  theme(legend.position="none")

plot

ggsave("Downloads/SRS404442.chrom-depth.svg", plot = plot, dpi = 200, units = "cm", height = 22, width = 22)