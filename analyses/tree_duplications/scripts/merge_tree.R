library(RRphylo)
library(manipulate)
library(ape)
library(phytools)
library(ggtree)
library(tidyverse)
library(RColorBrewer)
library(ggnewscale)

#### Metadata ####
metadata <- read.delim("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/data/derived/metadata_fixed.csv", header=TRUE, sep=",")

sublineage <- metadata %>%
                select(strain, vni_subdivision)%>%
                column_to_rownames("strain")
lineage <- metadata %>%
            select(strain, lineage)%>%
            column_to_rownames("strain")
dataset <- metadata %>%
            select(strain, dataset)%>%
            column_to_rownames("strain")
source <- metadata %>%
            select(strain, source)%>%
            column_to_rownames("strain")

#### Desjardins tree ####

desj_tree_path <- "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/data/raw/CryptoDiversity_Desjardins_Tree.tre"
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
write.tree(desj_tree, file = "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/data/processed/desj_tree.newick")

# Plot of Desjardins tree

d <- ggtree(desj_tree, layout = "circular") +  
            geom_tiplab(aes(label = label), size = 1.5, align =TRUE, linetype = "dashed", linesize = .05)+
            geom_treescale(x=0.1, y=0, width=0.1, offset = 0.5)
d1 <- gheatmap(d, lineage, width=.06, colnames=FALSE, offset=0.1) +
            scale_fill_brewer(palette = "Dark2", name="Lineage",  na.translate = FALSE)+
            new_scale_fill()
d2 <- gheatmap(d1, sublineage, width=.08, colnames=FALSE, offset=0.12) +
            scale_fill_brewer(palette = "Paired", name="SubLineage",  na.translate = FALSE)

ggsave("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/figures/tree_desjardins.png", d2, height = 10, width = 10, units = "in", dpi = 300)



#### Ashton tree ####

ashton_tree_path <- "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/data/raw/2017.06.09.all_ours_and_desj.snp_sites.mod.fa.cln.tree"
ashton_tree_unrooted <- read.tree(ashton_tree_path)

# Rename tips to use strains
ashton_tree_unrooted$tip.label <- sapply(ashton_tree_unrooted$tip.label, function(x) {
    if (x %in% metadata$run) {
        metadata$strain[metadata$run == x]
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


# Plot unrooted Ashton tree
pu <- ggtree(ashton_tree_unrooted, layout = "circular") +  
      geom_tiplab(aes(label = label), size = 0.5, align =TRUE, linetype = "dashed", linesize = .05)
pu1 <- gheatmap(pu, sublineage, width=.08, colnames=FALSE, offset=.01) +
    scale_fill_brewer(palette = "Paired", name="SubLineage",  na.translate = FALSE)

ggsave("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/figures/tree_ashton_unrooted.png", pu1, height = 10, width = 10, units = "in", dpi = 600)

# Root Ashton tree at the middle of the branch leading to VNIa
VNIa_root <- getMRCA(ashton_tree_unrooted, c("AD3-95a","Tu259-1"))
edge_length <- subset(ashton_tree_unrooted$edge.length, ashton_tree_unrooted$edge[,2] == VNIa_root)
ashton_tree <- reroot(ashton_tree_unrooted, VNIa_root, edge_length/2)

p <- ggtree(ashton_tree, layout = "circular") +  
      geom_tiplab(aes(label = label), size = 0.5, align =TRUE,linetype = "dashed", linesize = .05)
p1 <- gheatmap(p, sublineage, width=.08, colnames=FALSE, offset=.01) +
    scale_fill_brewer(palette = "Paired", name="SubLineage",  na.translate = FALSE)

ggsave("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/figures/tree_ashton.png", p1, height = 10, width = 10, units = "in", dpi = 600)
write.tree(ashton_tree, file = "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/data/processed/ashton_tree.newick")

# REctangular

ashton_tree <- read.tree("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/data/processed/ashton_tree.newick")
p <- ggtree(ashton_tree, layout = "rectangular") +  
      geom_tiplab(aes(label = label), size = 0.5, align =TRUE,linetype = "dashed", linesize = .05)
p1 <- gheatmap(p, sublineage, width=.08, colnames=FALSE, offset=.01) +
    scale_fill_brewer(palette = "Paired", name="SubLineage",  na.translate = FALSE)

ggsave("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/figures/tree_ashton_rectangular.pdf", p1, height = 10, width = 10, units = "in", dpi = 600)
#### Merge Desjardins and Ashton trees ####

# Specify clades in Desjardins tree

VNI <- c("Bt92", "Bt79")
VNI_node <- getMRCA(desj_tree, VNI)
VNII <- c("C2","C12")
VNII_node <- getMRCA(desj_tree, VNII)
VNB <- c("Bt7", "Bt34")
VNB_node <- getMRCA(desj_tree, VNB)

# Get the ages of the nodes from the original Desjardins tree
edge_lengths <- node.depth.edgelength(desj_tree)
node_labels <- c(desj_tree$tip.label, desj_tree$node.label)
edge_length_mapping <- data.frame(node = node_labels, edge_length = edge_lengths, max_length = max(edge_lengths))
edge_length_mapping <- edge_length_mapping %>% 
                        mutate(age = max_length - edge_length) %>%
                        rownames_to_column("node_id")
clade_ages <- edge_length_mapping %>% 
                filter(node_id %in% c(VNI_node, VNII_node, VNB_node))
nodeages <- c("Bt92-Bt79" = clade_ages$age[clade_ages$node_id == VNI_node],
             "C2-C12" = clade_ages$age[clade_ages$node_id == VNII_node],
             "Bt7-Bt34" = clade_ages$age[clade_ages$node_id == VNB_node])
tip_ages <- edge_length_mapping %>% 
                filter(node %in% metadata$strain)
tipages <- tip_ages$age
names(tipages) <- tip_ages$node

#Remove VNI clade from Desjardins tree to use it as backtree
VNI_tips <- tips(desj_tree, VNI_node)
backtree <- drop.tip(desj_tree, VNI_tips)
# Create the reference tables
reference <- data.frame(bind=c("CNS_289-BK8"),
                   reference=c("Bt7-Bt34"), # "H99"
                   poly=c(FALSE))

# Merge
merged <- tree.merger(backbone = backtree,
                        data=reference,
                        source.tree = ashton_tree,
                        plot=FALSE,
                        node.ages = nodeages,
                        tip.ages = tipages)
write.tree(merged, file = "/FastData/czirion/Crypto_Diversity_Pipeline/analyses/data/processed/merged_tree.newick")

# Plot merged tree
p <- ggtree(merged, layout = "circular", size = 0.1) +  
      geom_tiplab(aes(label = label), size = 0.5, align =TRUE, 
                    linetype = "dashed", linesize = .05)+
    geom_treescale(x=0.4, y=0, width=0.01, offset = 5)

m1 <- gheatmap(p, lineage, width=.05, colnames=FALSE, offset=.025) +
    scale_fill_brewer(palette = "Dark2", name="Lineage",  na.translate = FALSE)+
    new_scale_fill()

m2 <- gheatmap(m1, sublineage, width=.05, colnames=FALSE, offset=.042) +
    scale_fill_brewer(palette = "Paired", name="SubLineage",  na.translate = FALSE)

ggsave("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/figures/tree_merged_branchlengths.png", m2, height = 10, width = 10, units = "in", dpi = 600)

p <- ggtree(merged, layout = "circular", size = 0.1, branch.length = "none") +  
      geom_tiplab(aes(label = label), size = 0.5, align =TRUE, 
                    linetype = "dashed", linesize = .05)

m1 <- gheatmap(p, lineage, width=.05, colnames=FALSE, offset=3) +
    scale_fill_brewer(palette = "Dark2", name="Lineage",  na.translate = FALSE)+
    new_scale_fill()

m2 <- gheatmap(m1, sublineage, width=.05, colnames=FALSE, offset=5) +
    scale_fill_brewer(palette = "Paired", name="SubLineage",  na.translate = FALSE)

ggsave("/FastData/czirion/Crypto_Diversity_Pipeline/analyses/tree_duplications/results/figures/tree_merged.png", m2, height = 10, width = 10, units = "in", dpi = 600)



#### Explore the branch lengths of all trees and compare ####

## Merged tree

# Specify clades in merged tree
merged_VNI_node <- getMRCA(merged, VNI)
merged_VNII_node <- getMRCA(merged, VNII)
merged_VNB_node <- getMRCA(merged, VNB)
# Get the ages of the nodes from the original merged tree
merged_edge_lengths <- node.depth.edgelength(merged)
merged_node_labels <- c(merged$tip.label, merged$node.label)
merged_edge_length_mapping <- data.frame(node = merged_node_labels, edge_length = merged_edge_lengths, max_length = max(merged_edge_lengths))
merged_edge_length_mapping <- merged_edge_length_mapping %>% 
                                mutate(age = max_length - edge_length) %>%
                                rownames_to_column("node_id")
merged_clade_ages <- merged_edge_length_mapping %>% 
                        filter(node_id %in% c(merged_VNI_node, merged_VNII_node, merged_VNB_node))
merged_nodeages <- c("Bt92-Bt79" = merged_clade_ages$age[merged_clade_ages$node_id == merged_VNI_node],
                    "C2-C12" = merged_clade_ages$age[merged_clade_ages$node_id == merged_VNII_node],
                    "Bt7-Bt34" = merged_clade_ages$age[merged_clade_ages$node_id == merged_VNB_node])
tip_ages <- merged_edge_length_mapping %>% 
                filter(node %in% desj_tree$tip.label)
merged_tipages <- tip_ages$age
names(merged_tipages) <- tip_ages$node

merged_VNI_clade_length <- merged_edge_length_mapping %>% 
                    filter(node_id %in% merged_VNI_node)
merged_VNI_lengths <- merged_edge_length_mapping %>% 
                filter(node %in% VNI_tips)%>%
                mutate(VNI_length = VNI_clade_length$edge_length,
                        edge_length_VNI = edge_length - VNI_length)%>%
                select(node, edge_length_VNI)


## Desjardins tree

VNI_clade_length <- edge_length_mapping %>% 
                    filter(node_id %in% VNI_node)
VNI_lengths <- edge_length_mapping %>%
                filter(node %in% VNI_tips)%>%
                mutate(VNI_length = VNI_clade_length$edge_length,
                        edge_length_VNI = edge_length - VNI_length)%>%
                select(node, edge_length_VNI)

## Ashton ages

# Specify clades in Ashton tree
VNIa <- c("BK290", "CNS_289")
VNIa_node <- getMRCA(ashton_tree, VNIa)
VNIb <- c("C23", "AD3-41a")
VNIb_node <- getMRCA(ashton_tree, VNIb)
VNIc <- c("Bt11", "Bt20")
VNIc_node <- getMRCA(ashton_tree, VNIc)
VNI <- c("BK290", "Bt11")
VNI_node <- getMRCA(ashton_tree, VNI)
# Get the ages of the nodes from the original Ashton tree
ashton_tree_edge_lengths <- node.depth.edgelength(ashton_tree)
ashton_tree_node_labels <- c(ashton_tree$tip.label, ashton_tree$node.label)
ashton_tree_edge_length_mapping <- data.frame(node = ashton_tree_node_labels, edge_length = ashton_tree_edge_lengths, max_length = max(ashton_tree_edge_lengths))
ashton_tree_edge_length_mapping <- ashton_tree_edge_length_mapping %>% 
                                    mutate(age = max_length - edge_length) %>%
                                    rownames_to_column("node_id")
ashton_VNI_lengths <- ashton_tree_edge_length_mapping %>%
                        filter(node %in% VNI_tips)%>%
                        select(node, edge_length)

VNI_lengths_fuse <- merge(ashton_VNI_lengths,merged_VNI_lengths, by = "node")%>%
                    mutate(edge_length_diff = edge_length - edge_length_VNI)


ashton_tree_clade_ages <- ashton_tree_edge_length_mapping %>% 
                            filter(node_id %in% c(VNIa_node, VNIb_node, VNIc_node, VNI_node))
ashton_tree_nodeages <- c("BK290-CNS_289" = ashton_tree_clade_ages$age[ashton_tree_clade_ages$node_id == VNIa_node],
                        "C23-AD3-41a" = ashton_tree_clade_ages$age[ashton_tree_clade_ages$node_id == VNIb_node],
                        "Bt11-Bt20" = ashton_tree_clade_ages$age[ashton_tree_clade_ages$node_id == VNIc_node],
                        "BK290-Bt11" = ashton_tree_clade_ages$age[ashton_tree_clade_ages$node_id == VNI_node])

## Ashton unrooted ages

# Specify clades in Ashton tree
VNIa_node_unrooted <- getMRCA(ashton_tree_unrooted, VNIa)
VNIb_node_unrooted <- getMRCA(ashton_tree_unrooted, VNIb)
VNIc_node_unrooted <- getMRCA(ashton_tree_unrooted, VNIc)
VNI_node_unrooted <- getMRCA(ashton_tree_unrooted, VNI)
# Get the ages of the nodes from the original Ashton tree
unrooted_ashton_tree_edge_lengths <- node.depth.edgelength(ashton_tree_unrooted)
unrooted_ashton_tree_node_labels <- c(ashton_tree_unrooted$tip.label, ashton_tree_unrooted$node.label)
unrooted_ashton_tree_edge_length_mapping <- data.frame(node = unrooted_ashton_tree_node_labels, edge_length = unrooted_ashton_tree_edge_lengths, max_length = max(unrooted_ashton_tree_edge_lengths))
unrooted_ashton_tree_edge_length_mapping <- unrooted_ashton_tree_edge_length_mapping %>% 
                                                mutate(age = max_length - edge_length) %>%
                                                rownames_to_column("node_id")
unrooted_ashton_tree_clade_ages <- unrooted_ashton_tree_edge_length_mapping %>% 
                                    filter(node_id %in% c(VNIa_node_unrooted, VNIb_node_unrooted, VNIc_node_unrooted, VNI_node_unrooted))
unrooted_ashton_tree_nodeages <- c("BK290-CNS_289" = unrooted_ashton_tree_clade_ages$age[unrooted_ashton_tree_clade_ages$node_id == VNIa_node_unrooted],
                                "C23-AD3-41a" = unrooted_ashton_tree_clade_ages$age[unrooted_ashton_tree_clade_ages$node_id == VNIb_node_unrooted],
                                "Bt11-Bt20" = unrooted_ashton_tree_clade_ages$age[unrooted_ashton_tree_clade_ages$node_id == VNIc_node_unrooted],
                                "BK290-Bt11" = unrooted_ashton_tree_clade_ages$age[unrooted_ashton_tree_clade_ages$node_id == VNI_node_unrooted])
