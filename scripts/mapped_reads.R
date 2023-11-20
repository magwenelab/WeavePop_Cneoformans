log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)
metadata <- read.csv(snakemake@input[[2]], header = TRUE)%>%
    select(sample = Sample, strain = Strain, lineage = Group)

#metadata <- read.csv("./sample_metadata.csv", header = TRUE)%>%
#    select(sample = Sample, strain = Strain, lineage = Group)

stats <- read.delim(snakemake@input[[1]], sep =":", header = FALSE, col.names = c("stat", "value", "sample"))

stats <- read.delim("./read_stats.txt", sep =":", header = FALSE, col.names = c("stat", "value", "sample"))
stats <- stats %>% pivot_wider(names_from = stat, values_from = value)
colnames(stats) <- gsub(" ", "_", colnames(stats))
stats <- stats %>%
    mutate(total_reads = reads_mapped + reads_unmapped,
            percent_mapped = (reads_mapped/total_reads)*100)%>%
            as.data.frame()

stats <- left_join(stats, metadata, by="sample")

stats_long <- stats %>%
    select(sample, lineage, strain, mapped = percent_mapped, properly_paired = "percentage_of_properly_paired_reads_(%)") %>%
    pivot_longer(c(mapped, properly_paired), names_to = "measurement", values_to = "value")

plot <- ggplot(stats_long, aes(color = lineage, x=strain, y= value, shape = measurement))+
    geom_point()+
    ylim(0,100)+
    facet_wrap(~lineage,nrow = 1, scale = "free_x" )+
    scale_shape_discrete(labels=c('Mapped', 'Properly paired'))+
    labs(shape = NULL)+
    guides(color = FALSE)+ 
    theme_light()+
    xlab("Strain") +
    ylab("Percentage of reads") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#ggsave("./mapped_reads.svg", plot = plot, dpi = 200, units = "cm", height = 15, width = 15)
ggsave(snakemake@output[[1]], plot = plot, dpi = 200, units = "cm", height = 15, width = 15)