library(tidyverse)
library(RColorBrewer)
library(colorspace)


# This script creates named vectors of colors for different metadata fields.
# The objects here can be integrated in other R scripts with: source("metadata_colors.R")
# It can be used in ggplot plots with the scale_color_manual() and scale_fill_manual() functions.
# For example:
# scale_fill_manual(values = lineage_colors, name="Lineage")
# Assuming the fill aesthetic is used with the values in lineage_names.

# Chromosomes
chrom_colors <- c(colorRampPalette(brewer.pal(12, "Paired"))(14))
names(chrom_colors) <- c("chr01", "chr02", "chr03","chr04","chr05",
                         "chr06", "chr07", "chr08","chr09", "chr10",
                         "chr11", "chr12","chr13", "chr14")

# VNI sublineages
sublineage_names <- c("VNIb", "VNIc",
                    "VNIa-5", "VNIa-X",
                    "VNIa-32","VNIa-Y", 
                     "VNIa-93","VNIa-4")    
sublineage_colors <- c(brewer.pal(12, "Set3")[1:length(sublineage_names)])
names(sublineage_colors) <- sublineage_names     
sublineage_shading <- c("gray90", "gray70", 
                        "gray90", "gray70", 
                        "gray90", "gray70", 
                        "gray90", "gray70") 
names(sublineage_shading) <- sublineage_names  

# Lineages
lineage_names <- c("VNBI", "VNBII", "VNI", "VNII")
lineage_colors <- brewer.pal(8, "Dark2")[1:length(lineage_names)]
names(lineage_colors) <- lineage_names

# Datasets
dataset_names <- c("Ashton", "Desjardins", "Reference")
dataset_colors <- c(brewer.pal(3, "BuPu")[c(1, 3)], "white")
names(dataset_colors) <- dataset_names

# Source
source_names <- c("Clinical", "Environmental")
source_shapes <- c(16,17,0)
names(source_shapes) <- source_names

source_names <- c("Clinical", "Environmental")
source_colors <- brewer.pal(11, "BrBG")[c(9, 3)]
names(source_colors) <- source_names

# Mating type
mat_names <- c("α","a", NA)
# mat_colors <- c("black", "gray50", "gray90")
mat_colors <- c(brewer.pal(9, "Greens")[c(9,3)], "gray90")
names(mat_colors) <- mat_names

mat_shapes <- c(16,17,0)
names(mat_shapes) <- mat_names

# Countries
continent <- list(
    "North_America" = c("USA"),
    "South_America" = c("Argentina", "Brazil"),
    "Africa1" = c("Botswana", "Malawi", "S. Africa"), 
    "Africa2" = c("Tanzania", "Togo", "Uganda"),
    "Asia1" = c("China", "India", "Japan"), 
    "Asia2" = c("Laos", "Thailand", "Vietnam"),
    "Oceania" = c("Australia"),
    "Europe" = c("France", "Italy")
)

set2 <- brewer.pal(8, "Set2")
set2_dark <- darken(set2, 0.3)

country_colors <- unlist(Map(function(region, color, dark_color) {
    colorRampPalette(c(color, dark_color))(length(region))
}, continent, set2, set2_dark))

names(country_colors) <- unlist(continent, use.names = FALSE)

## Display colors
display_colors <- function(color_vector, title) {
    barplot(rep(1, length(color_vector)), col = color_vector, border = NA, 
                    names.arg = names(color_vector), las = 2, main = title)
}

# # Set up the plotting area
# par(mfrow = c(3, 2), mar = c(5, 4, 2, 1))
# # Display colors for each metadata field
# display_colors(dataset_colors, "Dataset Colors")
# display_colors(lineage_colors, "Lineage Colors")
# display_colors(sublineage_colors, "Sublineage Colors")
# display_colors(source_colors, "Source Colors")
# display_colors(mat_colors, "Mating Type Colors")
# display_colors(country_colors, "Country Colors")
# display_colors(chrom_colors, "Country Colors")

# #Reset plotting area
# par(mfrow = c(1, 1), mar = c(5, 4, 2, 1))

# display.brewer.all(colorblindFriendly = TRUE)
