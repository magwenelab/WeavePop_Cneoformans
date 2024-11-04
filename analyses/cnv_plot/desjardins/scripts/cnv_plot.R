suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ggtreeExtra))
suppressPackageStartupMessages(library(ape))
setwd("/FastData/czirion/Crypto_Desjardins/cnv_plot")

repeat_threshold <- 0.3
CHROM_NAMES_FILE <- "/FastData/czirion/Crypto_Desjardins/config/chromosomes.csv"
CHROM_LENGTHS_FILE <- "/FastData/czirion/Crypto_Desjardins/cnv_plot/chromosome_lengths.tsv"
LOCI_FILE <- "/FastData/czirion/Crypto_Desjardins/results_2024-08-19/references/loci_to_plot.tsv"
VNI_REPEATS_FILE <- "/FastData/czirion/Crypto_Desjardins/results_2024-08-19/references/VNI/repeats/VNI_repeats.bed"
VNII_REPEATS_FILE <- "/FastData/czirion/Crypto_Desjardins/results_2024-08-19/references/VNII/repeats/VNII_repeats.bed"
VNBI_REPEATS_FILE <- "/FastData/czirion/Crypto_Desjardins/results_2024-08-19/references/VNBI/repeats/VNBI_repeats.bed"
VNBII_REPEATS_FILE <- "/FastData/czirion/Crypto_Desjardins/results_2024-08-19/references/VNBII/repeats/VNBII_repeats.bed"
METADATA_FILE <- "/FastData/czirion/Crypto_Desjardins/config/metadata.csv"
CNV_FILE <- "/FastData/czirion/Crypto_Desjardins/results_2024-08-19/dataset/cnv/cnv_calls.tsv"

chrom_names <- read.delim(CHROM_NAMES_FILE, sep = ",", header = TRUE, col.names = c("Lineage", "Accession", "Chromosome"), colClasses = c("factor", "factor", "integer"))
lengths <- read.delim(CHROM_LENGTHS_FILE, sep = "\t", header = FALSE, col.names = c("Accession", "Length"), stringsAsFactors=TRUE)
loci <- read.delim(LOCI_FILE, sep = "\t", header = TRUE, stringsAsFactors=TRUE)
VNI <- read.delim(VNI_REPEATS_FILE, sep = "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_Type"), stringsAsFactors=TRUE)
VNII <- read.delim(VNII_REPEATS_FILE, sep = "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_Type"), stringsAsFactors=TRUE)
VNBI <- read.delim(VNBI_REPEATS_FILE, sep = "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_Type"), stringsAsFactors=TRUE)
VNBII <- read.delim(VNBII_REPEATS_FILE, sep = "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_Type"), stringsAsFactors=TRUE)
metadata <- read.delim(METADATA_FILE, sep = ",", header = TRUE)
cnv <- read.delim(CNV_FILE, sep = "\t", header = TRUE, stringsAsFactors=TRUE)

# =============================================================================
TREE_FILE <- "/FastData/czirion/Crypto_Desjardins/cnv_plot/CryptoDiversity_Desjardins_Tree.tre" 

tree <- read.tree(TREE_FILE)
outgroup_clade <- metadata %>%
    filter(lineage == "VNII")%>%
    select(strain)

outgroup_clade <- outgroup_clade$strain
tree <- root(tree, outgroup = outgroup_clade)

plot_tree <- ggplot(tree, branch.length = 'none', aes(x, y)) + 
    geom_tree() + 
    geom_tiplab(aes(label = label), size = 2) +
    theme_tree()

strain_order <- plot_tree$data %>%
    filter(isTip)%>%
    arrange(y)%>%
    pull(label)%>%
    rev()

# =============================================================================

# Order the rows of the metadata to match the strain_order without turning the column into a factor
metadata <- metadata[match(strain_order, metadata$strain),]
metadata$strain <- factor(metadata$strain, levels = strain_order)
metadata$sample <- factor(metadata$sample, levels = metadata$sample)

metadata <- metadata %>% 
    select(Sample = sample, Strain = strain, Lineage = lineage)
# =============================================================================

lengths <- left_join(lengths, chrom_names, by = "Accession")
lengths <- lengths %>% select(Chromosome, Lineage, End = Length)%>%
    group_by(Lineage)%>%
    mutate(Start = 0, Strain = "Chromosomes", Feature = "Chromosomes")
# =============================================================================

loci <- loci %>% rename("Accession" = seq_id)
loci <- left_join(loci, chrom_names, by = "Accession")
loci <- loci %>%
    select(Chromosome, Start = start, End = end, Feature = Loci, Lineage)

centromeres <- loci %>% filter(Feature == "Centromeres")%>%
    group_by(Lineage, Chromosome)%>%
    reframe(End = max(End), Start = min(Start), Feature = "Centromeres", Strain = "Centromeres")

mat <- loci %>% filter(Feature == "MAT")%>%
    group_by(Lineage, Chromosome)%>%
    reframe(End = max(End), Start = min(Start), Feature = "MAT Locus", Strain = "MAT Locus")
# =============================================================================

repeats <- bind_rows(lst(VNI, VNII, VNBI, VNBII), .id = "Lineage")
repeats <- left_join(repeats, chrom_names, by = c("Accession", "Lineage"))
repeats <- repeats %>%
    select(Chromosome, Start, End, Lineage)%>%
    group_by(Lineage)%>%
    mutate(Strain= "Repeats", Feature = "Repeats")

# =============================================================================

cnv <- left_join(cnv, metadata, by = "Sample")%>%
    filter(Repeat_fraction < repeat_threshold)
cnv <- left_join(cnv, chrom_names, by = c("Accession", "Lineage"))
cnv <- cnv %>%
    select(Chromosome, Start, End, Feature = Structure, Sample, Strain, Lineage)

# =============================================================================

df <- bind_rows(lst(cnv, mat, centromeres, repeats, lengths), .id = "Sample_Reference")
df$Sample_Reference <- ifelse(df$Sample_Reference == "cnv", "Samples", "Reference")
df$Sample_Reference <- factor(df$Sample_Reference, levels = c("Samples", "Reference"))
df$Feature <- factor(df$Feature, levels = c("DELETION", "DUPLICATION", "MAT Locus", "Centromeres", "Repeats", "Chromosomes"))
df$Strain <- factor(df$Strain, levels = c("Chromosomes", "Repeats","Centromeres", "MAT Locus",  levels(metadata$Strain)))
# =============================================================================

s_colors <- c("#E6AB02", "#66A61E", "#D95F02" ,"#1B9E77", "#A6761D" ,"#7570B3" )

# =============================================================================
lineages <- unique(df$Lineage)

for (lineage in lineages){
    df_lineage <- df %>% filter(Lineage == lineage)
    num_samples <- length(unique(df_lineage$Sample))
    pheight <- 5 + num_samples/2
    print(paste("Plotting lineage", lineage, "with", num_samples, "samples and size", pheight, "cm"))
    plot <- ggplot() +
        geom_segment(data = df_lineage, aes(x=Start, xend = End, y = Strain, yend = Strain, color = Feature), linewidth= 3) +
        facet_grid(Sample_Reference~Chromosome, scales = "free", space = "free", switch = "y") +
        scale_x_continuous(labels = scales::label_bytes(), breaks = function(x) max(x)) +
        scale_y_discrete(position = "right") +
        scale_color_manual(values = s_colors) +
        theme(panel.background = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 0, hjust = 1),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank())+
        labs(title = paste("Copy-number variants of lineage", lineage))
    print("Saving plot")
    ggsave(paste("cnv_", lineage, ".png", sep = ""), plot, height = pheight, width = 32, units = "cm", limitsize = FALSE)
}
