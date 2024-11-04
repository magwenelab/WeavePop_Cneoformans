suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ggtreeExtra))
suppressPackageStartupMessages(library(ape))
setwd("/FastData/czirion/ashton/cnv_plot")

repeat_threshold <- 0.3
CHROM_NAMES_FILE <- "/FastData/czirion/ashton/config/chromosome_names.csv"
CHROM_LENGTHS_FILE <- "/FastData/czirion/ashton/cnv_plot/chromosome_lengths.tsv"
LOCI_FILE <- "/FastData/czirion/ashton/results/references/loci_to_plot.tsv"
VNI_REPEATS_FILE <- "/FastData/czirion/ashton/results/references/VNI/repeats/VNI_repeats.bed"

METADATA_FILE <- "/FastData/czirion/ashton/config/sample_metadata.csv"
CNV_FILE <- "/FastData/czirion/ashton/results/dataset/cnv/cnv_calls.tsv"

chrom_names <- read.delim(CHROM_NAMES_FILE, sep = ",", header = TRUE, col.names = c("Lineage", "Accession", "Chromosome"), colClasses = c("factor", "factor", "integer"))
lengths <- read.delim(CHROM_LENGTHS_FILE, sep = "\t", header = FALSE, col.names = c("Accession", "Length"), stringsAsFactors=TRUE)
loci <- read.delim(LOCI_FILE, sep = "\t", header = TRUE, stringsAsFactors=TRUE)
VNI <- read.delim(VNI_REPEATS_FILE, sep = "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_Type"), stringsAsFactors=TRUE)
metadata <- read.delim(METADATA_FILE, sep = ",", header = TRUE)
cnv <- read.delim(CNV_FILE, sep = "\t", header = TRUE, stringsAsFactors=TRUE)

# =============================================================================
# TREE_FILE <- "/FastData/czirion/Crypto_Desjardins/cnv_plot/CryptoDiversity_Desjardins_Tree.tre" 

# tree <- read.tree(TREE_FILE)
# outgroup_clade <- metadata %>%
#     filter(lineage == "VNII")%>%
#     select(strain)

# outgroup_clade <- outgroup_clade$strain
# tree <- root(tree, outgroup = outgroup_clade)

# plot_tree <- ggplot(tree, branch.length = 'none', aes(x, y)) + 
#     geom_tree() + 
#     geom_tiplab(aes(label = label), size = 2) +
#     theme_tree()

# strain_order <- plot_tree$data %>%
#     filter(isTip)%>%
#     arrange(y)%>%
#     pull(label)%>%
#     rev()

# =============================================================================

# Order the rows of the metadata to match the strain_order without turning the column into a factor
metadata$strain <- as.factor(metadata$strain)
metadata$sample <- as.factor(metadata$sample)

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
    reframe(End = max(End), Start = min(Start), Feature = "MAT", Strain = "MAT")
# =============================================================================

repeats <- VNI
repeats$Lineage <- "VNI"
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

df <- bind_rows(lst(cnv, mat, centromeres, repeats, lengths), .id = "Things")
df$Feature <- factor(df$Feature, levels = c("DELETION", "DUPLICATION", "MAT", "Centromeres", "Repeats", "Chromosomes"))
df$Strain <- factor(df$Strain, levels = c("Chromosomes", "Repeats","Centromeres", "MAT",  levels(metadata$Strain)))
# =============================================================================

s_colors <- c("#E6AB02", "#66A61E", "#D95F02" ,"#1B9E77", "#A6761D" ,"#7570B3" )

# =============================================================================

chrom_plot <- ggplot() +
    geom_segment(data = df, aes(x=Start, xend = End, y = Strain, yend = Strain, color = Feature), linewidth= 3) +
    facet_grid(~Chromosome, scales = "free_x", space = "free_x") +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_color_manual(values = s_colors) +
    theme(panel.background = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_blank())+
    labs(title = "Copy-number variants of the Ashton dataset samples")
chrom_plot
ggsave("cnv.svg", chrom_plot, height = 150, width = 30, units = "cm", limitsize = FALSE)
