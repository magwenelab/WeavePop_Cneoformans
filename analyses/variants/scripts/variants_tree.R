library(ggtree)
library(TDbook)

## load `tree_nwk`, `df_info`, `df_alleles`, and `df_bar_data` from 'TDbook'
tree <- tree_nwk
snps <- df_alleles
snps_strainCols <- snps[1,] 
snps<-snps[-1,] # drop strain names
colnames(snps) <- snps_strainCols

gapChar <- "?"
snp <- t(snps)
lsnp <- apply(snp, 1, function(x) {
        x != snp[1,] & x != gapChar & snp[1,] != gapChar
    })
lsnp <- as.data.frame(lsnp)
lsnp$pos <- as.numeric(rownames(lsnp))
lsnp <- tidyr::gather(lsnp, name, value, -pos)
snp_data <- lsnp[lsnp$value, c("name", "pos")]

## visualize the tree 
p <- ggtree(tree) 

## attach the sampling information data set 
## and add symbols colored by location
p <- p %<+% df_info + geom_tippoint(aes(color=location))

## visualize SNP and Trait data using dot and bar charts,
## and align them based on tree structure
p + geom_facet(panel = "Trait", data = df_bar_data, geom = geom_col, 
                aes(x = dummy_bar_value, color = location, 
                fill = location), orientation = 'y', width = .6)

bar <- df_bar_data
