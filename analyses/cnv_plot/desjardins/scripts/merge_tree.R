library(RRphylo)
library(manipulate)
library(ape)
library(phytools)
library(ggtree)
library(tidyverse)
library(RColorBrewer)
library(ggnewscale)

#### Metadata ####
runs <- read.delim("/BigData/czirion/Crypto_Diversity_Pipeline/data/runs_table.tsv", header=TRUE, sep="\t")
metadata <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/results_joined_241204/02.Dataset/metadata.csv", header=TRUE, sep=",")
sublineage <- metadata %>%
    select(strain, vni_subdivision)%>%
    column_to_rownames("strain")
lineage <- metadata %>%
    select(strain, lineage)%>%
    column_to_rownames("strain")
#### Desjardins tree ####

desj_tree_path <- "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/data/CryptoDiversity_Desjardins_Tree.tre"
desj_tree <- read.tree(desj_tree_path)
# outgroup_clade <- metadata %>%
#     filter(lineage == "VNII")%>%
#     select(strain)

# outgroup_clade <- outgroup_clade$strain
# desj_tree <- root(desj_tree, outgroup = outgroup_clade)
# desj_tree <- root(desj_tree, node = 493)
# plot(desj_tree, show.tip.label = FALSE) 
# nodelabels()

# Reroot the tree at the middle of the branch leading to VNII
VNII_root <- getMRCA(desj_tree, c("C2","C12"))

edge_length <- subset(desj_tree$edge.length, desj_tree$edge[,2] == VNII_root)
desj_tree <- reroot(desj_tree, VNII_root, edge_length/2)


# Plot of Desjardins tree

d <- ggtree(desj_tree, layout = "circular") +  geom_tiplab(aes(label = label), size = 1.5, align =TRUE, 
                    linetype = "dashed", linesize = .05)
d1 <- gheatmap(d, lineage, width=.06, colnames=FALSE, offset=0.1) +
    scale_fill_brewer(palette = "Dark2", name="Lineage",  na.translate = FALSE)+
    new_scale_fill()
d2 <- gheatmap(d1, sublineage, width=.08, colnames=FALSE, offset=0.12) +
    scale_fill_brewer(palette = "Paired", name="SubLineage",  na.translate = FALSE)

ggsave("analyses/cnv_plot/desjardins/results/desjardins_tree.png", d2, height = 10, width = 10, units = "in", dpi = 300)



#### Ashton tree ####

ashton_tree_path <- "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/cnv_plot/desjardins/data/2017.06.09.all_ours_and_desj.snp_sites.mod.fa.cln.tree"

ashton_tree_unrooted <- read.tree(ashton_tree_path)
ashton_tree_unrooted$tip.label <- sapply(ashton_tree_unrooted$tip.label, function(x) {
    if (x %in% runs$run) {
        runs$sample[runs$run == x]
    } else {
        x
    }
})
ashton_tree_unrooted$tip.label <- sapply(ashton_tree_unrooted$tip.label, function(x) {
    if (x %in% metadata$sample) {
        metadata$strain[metadata$sample == x]
    } else {
        x
    }
})
ashton_tree_unrooted$tip.label <- sapply(ashton_tree_unrooted$tip.label, function(x) {
    if (x %in% metadata$name) {
        metadata$strain[metadata$name == x]
    } else {
        x
    }
})


# Root Ashton tree at the middle of the branch leading to VNIa
VNIa_root <- getMRCA(ashton_tree_unrooted, c("AD3-95a","Tu259-1"))

edge_length <- subset(ashton_tree_unrooted$edge.length, ashton_tree_unrooted$edge[,2] == VNIa_root)
ashton_tree <- reroot(ashton_tree_unrooted, VNIa_root, edge_length/2)

p <- ggtree(ashton_tree, layout = "circular") +  
      geom_tiplab(aes(label = label), size = 1, align =TRUE, 
                    linetype = "dashed", linesize = .05)
p1 <- gheatmap(p, sublineage, width=.08, colnames=FALSE, offset=.03) +
    scale_fill_brewer(palette = "Paired", name="SubLineage",  na.translate = FALSE)

ggsave("analyses/cnv_plot/desjardins/results/ashton_tree.png", p1, height = 10, width = 10, units = "in", dpi = 300)

## Specify clades in Desjardins tree

VNI <- c("Tu241-1","C23", "Bt43", "AD3-95a" )
VNI_node <- getMRCA(desj_tree, VNI)
VNII <- c("C2","C12")
VNII_node <- getMRCA(desj_tree, VNII)
VNB <- c("Bt7", "Bt34")
VNB_node <- getMRCA(desj_tree, VNB)

# Get the ages of the nodes
edge_lengths <- node.depth.edgelength(desj_tree)
node_labels <- c(desj_tree$tip.label, desj_tree$node.label)
edge_length_mapping <- data.frame(node = node_labels, edge_length = edge_lengths, max_length = max(edge_lengths))
edge_length_mapping <- edge_length_mapping %>% 
    mutate(age = max_length - edge_length) %>%
    rownames_to_column("node_id")

clade_ages <- edge_length_mapping %>% 
  filter(node_id %in% c(VNI_node, VNII_node, VNB_node))



#Remove VNI clade from Desjardins tree
VNI_root_des <- getMRCA(desj_tree, c("Tu241-1","C23", "Bt43", "AD3-95a" ) )

# VNI_tips <- tips(desj_tree, VNI_root_des)
# backtree <- drop.tip(desj_tree, c(VNI_tips[-which(tips(desj_tree,VNI_root_des)%in%
#                                              c("Tu241-1","C23", "Bt43", "AD3-95a" ))]))



# Define VNI clade from Ashton tree
# c("Tu241-1","VNIb", "VNIc", "VNIa" )
VNI_root <- getMRCA(ashton_tree, c("Tu241-1","C23", "Bt43", "AD3-95a" ))
ashton_tree$node.label[(VNI_root-Ntip(ashton_tree))]<-"VNI" # assigning node labels

VNI_root <- getMRCA(desj_tree, c("Tu241-1","C23", "Bt43", "AD3-95a" ))
# Remove all VNI tips except H99
VNI_tips <- tips(desj_tree, VNI_root)
backtree <- drop.tip(desj_tree, VNI_tips[-which(tips(desj_tree,VNI_root) %in% c("H99"))])
# Backtree Ashton
# Source Desjardins
# Create the reference table

reference <- data.frame(bind=c("CNS_289-BK8"), # Tips of VNI in Ashton tree, VNI-VNB tips in merged tree
                   reference=c("H99"), # Tips of VNBI and VNBII in Desjarins tree, VNII tips in Desjardins tree
                   poly=c(FALSE))

# Merge
merged <- tree.merger(backbone = desj_tree,data=reference,source.tree = ashton_tree,plot=FALSE)
plot(merged)
merged <- drop.tip(merged, tips(merged, "H99"))


p <- ggtree(merged, layout = "circular", size = 0.1) +  
      geom_tiplab(aes(label = label), size = .5, align =TRUE, 
                    linetype = "dashed", linesize = .05)

m1 <- gheatmap(p, sublineage, width=.08, colnames=FALSE, offset=.01) +
    scale_fill_brewer(palette = "Paired", name="SubLineage",  na.translate = FALSE)+
    new_scale_fill()

m2 <- gheatmap(m1, lineage, width=.06, colnames=FALSE, offset=.02) +
    scale_fill_brewer(palette = "Dark2", name="Lineage",  na.translate = FALSE)

ggsave("analyses/cnv_plot/desjardins/results/merged_tree.png", m2, height = 10, width = 10, units = "in", dpi = 300)


VNI_a <- getMRCA(merged, c("BK43", "BK224"))
merged_minimal <- drop.tip(merged, tips(merged, VNI_a))
ggtree(merged_minimal, layout="circular")
