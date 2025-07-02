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
vnisublineage_names <- c("VNIb", "VNIc",
                    "VNIa-5", "VNIa-X",
                    "VNIa-32","VNIa-Y", 
                     "VNIa-93","VNIa-4")    
sublineage_colors <- c(brewer.pal(12, "Set3")[1:length(vnisublineage_names)])
names(sublineage_colors) <- vnisublineage_names     
vnisublineage_shading <- c("gray90", "gray70", 
                        "gray90", "gray70", 
                        "gray90", "gray70", 
                        "gray90", "gray70") 
names(vnisublineage_shading) <- vnisublineage_names  
#  sublineages
sublineage_names <- c("VNII", "VNBII", "VNBI", 
                     "VNIb", "VNIc",
                    "VNIa-5", "VNIa-X",
                    "VNIa-32","VNIa-Y", 
                     "VNIa-93","VNIa-4"
                     )    
sublineage_shading <- c("gray70", "gray90", "gray70",
                        "gray90", "gray70", 
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
    "Australasia" = c("Australia"),
    "Europe" = c("France", "Italy")
)

set2 <- brewer.pal(8, "Set2")
set2_dark <- darken(set2, 0.3)

country_colors <- unlist(Map(function(region, color, dark_color) {
    colorRampPalette(c(color, dark_color))(length(region))
}, continent, set2, set2_dark))

names(country_colors) <- unlist(continent, use.names = FALSE)


# Contienents

continents <- c("North America", "South America",
                "Europe", "Africa", 
                "Asia", "Australasia")
continent_colors <- brewer.pal(6, "Dark2")
names(continent_colors) <- continents

# Category of duplications by size
# category_size_colors <-c("#8E0152","#2171B5", "#DEEBF7","#F7FBFF" )
# names(category_size_colors) <- c("Absent","Small", "Medium", "Large")

# # Category of duplications by aneuploidy
# category_aneuploidy_colors <-c("#8E0152","#2171B5", "#DEEBF7" )
# names(category_aneuploidy_colors) <- c("Full", "Partial", "Euploid")

## Display colors
display_colors <- function(color_vector, title) {
    barplot(rep(1, length(color_vector)), col = color_vector, border = NA, 
            names.arg = names(color_vector), las = 1, main = title, yaxt = "n")
}


# display_colors(mat_colors, "Mating Type Colors")

# # Set up the plotting area
# par(mfrow = c(3, 1), mar = c(3, 3, 3, 3))
# # Display colors for each metadata field
# display_colors(dataset_colors, "Dataset Colors")
# display_colors(lineage_colors, "Lineage Colors")
# display_colors(sublineage_colors, "Sublineage Colors")
# display_colors(source_colors, "Source")
# display_colors(category_aneuploidy_colors, "Aneuploidy")
# display_colors(continent_colors, "Continent")
# display_colors(country_colors, "Country Colors")
# display_colors(chrom_colors, "Country Colors")

#Reset plotting area
# par(mfrow = c(8, 2), mar = c(1, 1, 1, 1))

# Display colorblind-friendly palettes from RColorBrewer
# display.brewer.all(colorblindFriendly = TRUE)

# Plot colorblind-friendly palettes for 6 categorical values
# Plot colorblind-friendly palettes for 6 categorical values from palettes not available in RColorBrewer

# Use palettes from the colorspace package that are colorblind-friendly and not in RColorBrewer

# List of colorblind-friendly palettes from colorspace (not in RColorBrewer)
# cs_palettes <- c("Okabe-Ito", "HCL_purple", "HCL_cyan", "HCL_orange", "HCL_red", "HCL_green", "HCL_blue")

# par(mfrow = c(length(cs_palettes), 1), mar = c(2, 6, 2, 2))
# for (pal in cs_palettes) {
#     if (pal == "Okabe-Ito") {
#         colors <- qualitative_hcl(6, palette = "Okabe-Ito")
#     } else {
#         colors <- qualitative_hcl(6, palette = pal)
#     }
#     barplot(rep(1, 6), col = colors, border = NA, names.arg = colors, las = 1, 
#                     main = paste("colorspace Palette:", pal), yaxt = "n")
# }

# Get all palette names
# palettes <- brewer.pal.info

# # Function to print palette colors with hex codes
# for (pal in rownames(palettes)) {
#   n <- palettes[pal, "maxcolors"]
#   colors <- brewer.pal(n, pal)
#   cat(sprintf("\nPalette: %s\n", pal))
#   print(setNames(colors, colors))
# }
