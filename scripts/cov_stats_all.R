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

global <-read.csv("results/coverage_global_good.csv", header = FALSE, col.names = c("Global_Mean", "Global_Median", "Sample"), stringsAsFactors = TRUE)
chromosome <-read.csv("results/coverage_good.csv", header = FALSE, col.names = c("Chromosome", "Measurement", "Value", "Sample"), stringsAsFactors = TRUE)
global <- left_join(global, metadata, by = "Sample")
chromosome$Chromosome <- as.factor(chromosome$Chromosome)

chromosome <- pivot_wider(chromosome, names_from = Measurement, values_from = Value)
chromosome <- left_join(chromosome, global, by = "Sample") 
chromosome <- chromosome%>%
    group_by(Chromosome, Sample)%>%
    mutate(pmean= round(Mean/Global_Mean, 2))%>%
    mutate(pmedian= round(Median/Global_Median, 2))%>%
    ungroup()

write_csv(chromosome, "results/proprotional_coverage_good.csv", col_names = TRUE)


# Global
topylim <- max(global$Global_Mean) + max(global$Global_Mean/10)
color_stat = c("Mean" = "#008837", "Median" = "#7b3294")

g <- ggplot(global, aes(x=reorder(name, -Global_Mean, sum)))+
    geom_point(aes(y= Global_Mean, color = "Mean"))+
    geom_point(aes(y= Global_Median, color = "Median"))+
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

toplim <- ceiling(max(chromosome$pmedian))
values <- seq(0, toplim, by = 1)
colors <- brewer.pal(n = toplim +1, name = "Dark2") 
ylabel <- "Ploidy"

medianplot <- ggplot(chromosome, aes(x=reorder(name, -Global_Mean, sum), y= pmedian))+
    geom_point(aes(color= pmedian))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(Group))+
    scale_color_gradientn(colors = colors, breaks = values,limits = c(0, toplim), values = rescale(values), guide = "colorbar", name = ylabel) +
    theme_light()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
          panel.grid.minor = element_blank())+
    labs(title = "Proportional median coverage of chromosome",
         x = "Sample",
         y = ylabel)

ggsave("./results/cov_prop_median_good.svg", plot = medianplot, dpi = 50, units = "cm", height = 30, width = 60)

# Mean by Chromosome

toplim <- ceiling(max(chromosome$pmean))
values <- seq(0, toplim, by = 1)
colors <- brewer.pal(n = toplim +1, name = "Dark2") 
ylabel <- "Ploidy"

meanplot <- ggplot(chromosome, aes(x=reorder(name, -Global_Mean, sum), y= pmean))+
    geom_point(aes(color= pmean))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(Group))+
    scale_color_gradientn(colors = colors, breaks = values,limits = c(0, toplim), values = rescale(values), guide = "colorbar", name = ylabel) +
    theme_light()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
          panel.grid.minor = element_blank())+
    labs(title = "Proportional mean coverage of chromosome",
         x = "Sample",
         y = ylabel)

ggsave("./results/cov_prop_mean_good.svg", plot = meanplot, dpi = 50, units = "cm", height = 30, width = 60)

#### All mappings ####

global <-read.csv("results/coverage_global_raw.csv", header = FALSE, col.names = c("Global_Mean", "Global_Median", "Sample"), stringsAsFactors = TRUE)
chromosome <-read.csv("results/coverage_raw.csv", header = FALSE, col.names = c("Chromosome", "Measurement", "Value", "Sample"), stringsAsFactors = TRUE)
global <- left_join(global, metadata, by = "Sample")
chromosome$Chromosome <- as.factor(chromosome$Chromosome)

chromosome <- pivot_wider(chromosome, names_from = Measurement, values_from = Value)
chromosome <- left_join(chromosome, global, by = "Sample") 
chromosome <- chromosome%>%
    group_by(Chromosome, Sample)%>%
    mutate(pmean= round(Mean/Global_Mean, 2))%>%
    mutate(pmedian= round(Median/Global_Median, 2))%>%
    ungroup()

write_csv(chromosome, "results/proprotional_coverage_raw.csv", col_names = TRUE)

# Global
topylim <- max(global$Global_Mean) + max(global$Global_Mean/10)
color_stat = c("Mean" = "#008837", "Median" = "#7b3294")

g <- ggplot(global, aes(x=reorder(name, -Global_Mean, sum)))+
    geom_point(aes(y= Global_Mean, color = "Mean"))+
    geom_point(aes(y= Global_Median, color = "Median"))+
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

toplim <- ceiling(max(chromosome$pmedian))
values <- seq(0, toplim, by = 1)
colors <- brewer.pal(n = toplim +1, name = "Dark2") 
ylabel <- "Ploidy"

medianplot <- ggplot(chromosome, aes(x=reorder(name, -Global_Mean, sum), y= pmedian))+
    geom_point(aes(color= pmedian))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(Group))+
    scale_color_gradientn(colors = colors, breaks = values,limits = c(0, toplim), values = rescale(values), guide = "colorbar", name = ylabel) +
    theme_light()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
          panel.grid.minor = element_blank())+
    labs(title = "Proportional median coverage of chromosome",
         x = "Sample",
         y = ylabel)

ggsave("./results/cov_prop_median_raw.svg", plot = medianplot, dpi = 50, units = "cm", height = 30, width = 60)

# Mean by Chromosome

toplim <- ceiling(max(chromosome$pmean))
values <- seq(0, toplim, by = 1)
colors <- brewer.pal(n = toplim +1, name = "Dark2") 
ylabel <- "Ploidy"

meanplot <- ggplot(chromosome, aes(x=reorder(name, -Global_Mean, sum), y= pmean))+
    geom_point(aes(color= pmean))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(Group))+
    scale_color_gradientn(colors = colors, breaks = values,limits = c(0, toplim), values = rescale(values), guide = "colorbar", name = ylabel) +
    theme_light()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
          panel.grid.minor = element_blank())+
    labs(title = "Proportional mean coverage of chromosome",
         x = "Sample",
         y = ylabel)

ggsave("./results/cov_prop_mean_raw.svg", plot = meanplot, dpi = 50, units = "cm", height = 30, width = 60)