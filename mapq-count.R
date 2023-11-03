log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)

print("Reading TSV file")
mapq<- read.delim(snakemake@input[[1]], header = FALSE, col.names = c("MAPQ", "Count"))

print("Plotting MAPQ count")

plot <- ggplot(mapq, aes(x=MAPQ, y=Count))+
  geom_col(fill = "hotpink4")+
  scale_y_continuous(name = "Count", labels = comma)+
  theme_light()

print("Saving plot")
ggsave(snakemake@output[[1]], plot = plot, dpi = 200, units = "cm", height = 10, width = 10)
print("Done!")