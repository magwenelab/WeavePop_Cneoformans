log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)
library(ggnewscale)

global_raw <-read.csv("results/coverage_global_raw.csv", header = FALSE, col.names = c("Mean", "Median", "Sample"), stringsAsFactors = TRUE)
global_good <-read.csv("results/coverage_global_good.csv", header = FALSE, col.names = c("Mean", "Median", "Sample"), stringsAsFactors = TRUE)
raw <-read.csv("results/coverage_raw.csv", header = FALSE, col.names = c("Chromosome", "Measurement", "Value", "Sample"), stringsAsFactors = TRUE)
metadata <- read.csv("sample_metadata.csv", header = TRUE, stringsAsFactors = TRUE)

# GLOBAL RAW
global_raw <- left_join(global_raw, metadata, by = "Sample")
toplim <- max(global_raw$Mean) + max(global_raw$Mean/10)
plot <- ggplot(global_raw, aes(x=reorder(Sample, -Mean, sum)))+
    geom_point(aes(y= Mean), color = "blue")+
    ylim(0,toplim)+
    facet_grid(~Group,scale = "free_x" , space='free_x')+
    theme_light()+
    xlab("Sample") +
    ylab("Mean global coverage") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))

ggsave("./results/cov_raw_all.svg", plot = plot, dpi = 50, units = "cm", height = 15, width = 60)

# BY CHROMOSOME RAW

raw$Chromosome <- as.factor(raw$Chromosome)
means <- filter(raw, Measurement == "Mean") %>% select(-Measurement)
medians <- filter(raw, Measurement == "Median") %>% select(-Measurement)

means <- left_join(means, global_raw, by = "Sample") 
means <- means%>%
    group_by(Chromosome, Sample)%>%
    mutate(pmean = Value/Mean)
colors <- c("blue", "green", "yellow", "red")
toplim <- max(means$pmean) + max(means$pmean/10)
plot <- ggplot(means, aes(x=Sample, y= pmean))+
    geom_point(aes(color= pmean))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(Group))+
    scale_color_gradientn(colors = colors, values = rescale(c(0, 1, 2, 3)), guide = "none") +
    guides(color = guide_colorbar(barwidth = 1, barheight = 4, nbin = 100, title.position = "top")) +
    #scale_color_discrete(labels=c('Mean', 'Median'),guide = guide_legend(title = NULL))+
    #labs(shape = NULL)+
    theme_light()+
    xlab("Sample") +
    ylab("Ploidy (proportional mean coverage of chromosome)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),panel.grid.minor = element_blank(),)

ggsave("./results/cov_prop_mean.svg", plot = plot, dpi = 50, units = "cm", height = 30, width = 60)


# GLOBAL GOOD
global_good <- left_join(global_good, metadata, by = "Sample")
toplim <- max(global_good$Mean) + max(global_good$Mean/10)
g <- ggplot(global_good, aes(x=reorder(Sample, -Mean, sum)))+
    geom_point(aes(y= Mean), color = "darkgreen")+
    ylim(0,toplim)+
    facet_grid(~Group,scale = "free_x" , space='free_x')+
    theme_light()+
    xlab("Sample") +
    ylab("Mean global coverage") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))

ggsave("./results/cov_good_all.svg", plot = g, dpi = 50, units = "cm", height = 30, width = 60)