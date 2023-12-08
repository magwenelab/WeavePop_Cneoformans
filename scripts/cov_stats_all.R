log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)
library(ggnewscale)

global_good <-read.csv("results/coverage_global_good.csv", header = FALSE, col.names = c("Mean", "Median", "Sample"), stringsAsFactors = TRUE)
global_raw <-read.csv("results/coverage_global_raw.csv", header = FALSE, col.names = c("Mean", "Median", "Sample"), stringsAsFactors = TRUE)
good <-read.csv("results/coverage_good.csv", header = FALSE, col.names = c("Chromosome", "Measurement", "Value", "Sample"), stringsAsFactors = TRUE)
raw <-read.csv("results/coverage_raw.csv", header = FALSE, col.names = c("Chromosome", "Measurement", "Value", "Sample"), stringsAsFactors = TRUE)

metadata <- read.csv("sample_metadata.csv", header = TRUE, stringsAsFactors = TRUE)
metadata <- mutate(metadata, name = paste(Strain, Sample, sep=" " ))

global_good <- left_join(global_good, metadata, by = "Sample")
global_raw <- left_join(global_raw, metadata, by = "Sample")

good$Chromosome <- as.factor(good$Chromosome)
raw$Chromosome <- as.factor(raw$Chromosome)


means <- filter(good, Measurement == "Mean") %>% select(-Measurement)
means <- filter(good, Measurement == "Median") %>% select(-Measurement)
means <- left_join(means, global_good, by = "Sample") 
means <- means%>%
    group_by(Chromosome, Sample)%>%
    mutate(proportion = Value/Mean)

medians <- filter(good, Measurement == "Median") %>% select(-Measurement)
medians <- left_join(medians, global_good, by = "Sample") 
medians <- medians%>%
    group_by(Chromosome, Sample)%>%
    mutate(proportion = Value/Median)



# GLOBAL GOOD
toplim <- max(global_good$Mean) + max(global_good$Mean/10)

g <- ggplot(global_good, aes(x=reorder(name, -Mean, sum)))+
    geom_point(aes(y= Mean), color = "aquamarine4")+
    ylim(0,toplim)+
    facet_grid(~Group,scale = "free_x" , space='free_x')+
    theme_light()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))+
    labs(title= "Mean global coverage of good mapping quality reads",
            x= "Sample",
            y= "Coverage (X)")

ggsave("./results/cov_good_all.svg", plot = g, dpi = 50, units = "cm", height = 30, width = 60)

# Good Median by Chromosome


toplim <- ceiling(max(medians$proportion))
values <- seq(0, toplim, by = 1)
colors <- brewer.pal(n = toplim +1, name = "Dark2") 
ylabel <- "Ploidy"

plot <- ggplot(medians, aes(x=reorder(name, -Mean, sum), y= proportion))+
    geom_point(aes(color= proportion))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(Group))+
    scale_color_gradientn(colors = colors, breaks = values,limits = c(0, ceiling(max(medians$proportion))), values = rescale(values), guide = "colorbar", name = ylabel) +
    theme_light()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
          panel.grid.minor = element_blank())+
    labs(title = "Proportional median coverage of chromosome",
         x = "Sample",
         y = ylabel)

ggsave("./results/cov_prop_median_good.svg", plot = plot, dpi = 50, units = "cm", height = 30, width = 60)