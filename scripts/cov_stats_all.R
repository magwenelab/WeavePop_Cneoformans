#log <- file(snakemake@log[[1]], open="wt")
#sink(log, type = "output")
#sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)
library(ggnewscale)

metadata <- read.csv("sample_metadata.csv", header = TRUE, stringsAsFactors = TRUE)
metadata <- mutate(metadata, name = paste(Strain, Sample, sep=" " ))

#### Good quality mappings ####

global <-read.csv("results/coverage_global_good.csv", header = FALSE, col.names = c("Mean", "Median", "Sample"), stringsAsFactors = TRUE)
chromosome <-read.csv("results/coverage_good.csv", header = FALSE, col.names = c("Chromosome", "Measurement", "Value", "Sample"), stringsAsFactors = TRUE)
global <- left_join(global, metadata, by = "Sample")
chromosome$Chromosome <- as.factor(chromosome$Chromosome)

means <- filter(chromosome, Measurement == "Mean") %>% select(-Measurement)
means <- left_join(means, global, by = "Sample") 
means <- means%>%
    group_by(Chromosome, Sample)%>%
    mutate(proportion = Value/Mean)

medians <- filter(chromosome, Measurement == "Median") %>% select(-Measurement)
medians <- left_join(medians, global, by = "Sample") 
medians <- medians%>%
    group_by(Chromosome, Sample)%>%
    mutate(proportion = Value/Median)

# Global
topylim <- max(global$Mean) + max(global$Mean/10)
color_stat = c("Mean" = "#008837", "Median" = "#7b3294")

g <- ggplot(global, aes(x=reorder(name, -Mean, sum)))+
    geom_point(aes(y= Mean, color = "Mean"))+
    geom_point(aes(y= Median, color = "Median"))+
    scale_color_manual(values= color_stat, name = "")+ 
    ylim(0,topylim)+
    facet_grid(~Group,scale = "free_x" , space='free_x')+
    theme_light()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))+
    labs(title= "Global coverage",
            x= "Sample",
            y= "Coverage (X)")

ggsave("./results/cov_good_all.svg", plot = g, dpi = 50, units = "cm", height = 30, width = 60)

# Median by Chromosome

toplim <- ceiling(max(medians$proportion))
values <- seq(0, toplim, by = 1)
colors <- brewer.pal(n = toplim +1, name = "Dark2") 
ylabel <- "Ploidy"

medianplot <- ggplot(medians, aes(x=reorder(name, -Mean, sum), y= proportion))+
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

ggsave("./results/cov_prop_median_good.svg", plot = medianplot, dpi = 50, units = "cm", height = 30, width = 60)

# Mean by Chromosome

toplim <- ceiling(max(means$proportion))
values <- seq(0, toplim, by = 1)
colors <- brewer.pal(n = toplim +1, name = "Dark2") 
ylabel <- "Ploidy"

meanplot <- ggplot(means, aes(x=reorder(name, -Mean, sum), y= proportion))+
    geom_point(aes(color= proportion))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(Group))+
    scale_color_gradientn(colors = colors, breaks = values,limits = c(0, ceiling(max(means$proportion))), values = rescale(values), guide = "colorbar", name = ylabel) +
    theme_light()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
          panel.grid.minor = element_blank())+
    labs(title = "Proportional mean coverage of chromosome",
         x = "Sample",
         y = ylabel)

ggsave("./results/cov_prop_mean_good.svg", plot = meanplot, dpi = 50, units = "cm", height = 30, width = 60)

#### All mappings ####

global <-read.csv("results/coverage_global_raw.csv", header = FALSE, col.names = c("Mean", "Median", "Sample"), stringsAsFactors = TRUE)
chromosome <-read.csv("results/coverage_raw.csv", header = FALSE, col.names = c("Chromosome", "Measurement", "Value", "Sample"), stringsAsFactors = TRUE)
global <- left_join(global, metadata, by = "Sample")
chromosome$Chromosome <- as.factor(chromosome$Chromosome)

means <- filter(chromosome, Measurement == "Mean") %>% select(-Measurement)
means <- left_join(means, global, by = "Sample") 
means <- means%>%
    group_by(Chromosome, Sample)%>%
    mutate(proportion = Value/Mean)

medians <- filter(chromosome, Measurement == "Median") %>% select(-Measurement)
medians <- left_join(medians, global, by = "Sample") 
medians <- medians%>%
    group_by(Chromosome, Sample)%>%
    mutate(proportion = Value/Median)

# Global
topylim <- max(global$Mean) + max(global$Mean/10)
color_stat = c("Mean" = "#008837", "Median" = "#7b3294")

g <- ggplot(global, aes(x=reorder(name, -Mean, sum)))+
    geom_point(aes(y= Mean, color = "Mean"))+
    geom_point(aes(y= Median, color = "Median"))+
    scale_color_manual(values= color_stat, name = "")+ 
    ylim(0,topylim)+
    facet_grid(~Group,scale = "free_x" , space='free_x')+
    theme_light()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))+
    labs(title= "Global coverage",
            x= "Sample",
            y= "Coverage (X)")

ggsave("./results/cov_raw_all.svg", plot = g, dpi = 50, units = "cm", height = 30, width = 60)

# Median by Chromosome

toplim <- ceiling(max(medians$proportion))
values <- seq(0, toplim, by = 1)
colors <- brewer.pal(n = toplim +1, name = "Dark2") 
ylabel <- "Ploidy"

medianplot <- ggplot(medians, aes(x=reorder(name, -Mean, sum), y= proportion))+
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

ggsave("./results/cov_prop_median_raw.svg", plot = medianplot, dpi = 50, units = "cm", height = 30, width = 60)

# Mean by Chromosome

toplim <- ceiling(max(means$proportion))
values <- seq(0, toplim, by = 1)
colors <- brewer.pal(n = toplim +1, name = "Dark2") 
ylabel <- "Ploidy"

meanplot <- ggplot(means, aes(x=reorder(name, -Mean, sum), y= proportion))+
    geom_point(aes(color= proportion))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(Group))+
    scale_color_gradientn(colors = colors, breaks = values,limits = c(0, ceiling(max(means$proportion))), values = rescale(values), guide = "colorbar", name = ylabel) +
    theme_light()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
          panel.grid.minor = element_blank())+
    labs(title = "Proportional mean coverage of chromosome",
         x = "Sample",
         y = ylabel)

ggsave("./results/cov_prop_mean_raw.svg", plot = meanplot, dpi = 50, units = "cm", height = 30, width = 60)