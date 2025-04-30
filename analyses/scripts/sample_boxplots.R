windows_path <- "../Crypto_Ashton/results/04.Intermediate_files/01.Samples/depth_quality/ERS1142845/depth_by_windows.tsv"
chromosomes_path <- "../Crypto_Ashton/results/04.Intermediate_files/03.References/chromosome_lengths.tsv"
windows <- read.delim(
                windows_path,
                header = FALSE,
                col.names = c("accession", "start", "end", "depth", "norm_depth", "smooth_depth"),
                sep = "\t",
                stringsAsFactors = TRUE)
chromosomes <- read.delim(
                    chromosomes_path,
                    header = TRUE,
                    sep = ",",
                    stringsAsFactors = TRUE) %>%
                mutate(chromosome = paste("chr", str_pad(chromosome, 2, pad = "0"), set = ""))%>%
                filter(accession %in% unique(windows$accession))%>%
                arrange(length)


chromosomes$chromosome <- factor(chromosomes$chromosome, levels = unique(chromosomes$chromosome))

windows <- left_join(windows, chromosomes, by = "accession")

c <- ggplot(windows, aes(x = chromosome, y = norm_depth))+
    geom_quasirandom(alpha = 0.05)+
    geom_boxplot(aes(color = chromosome), outlier.shape = NA)+
    scale_y_continuous(limits = c(0,4), breaks = seq(0,5, by = 0.5))+
    theme_minimal()+
    theme(legend.position = "none")+
    labs(y="Normalized Depth of Windows",
        x = "Chromosomes Ordered by Length",
        title = "Distribution of Normalized, Depth of Windows per Chromosome")
g <- ggplot(windows, aes(y = norm_depth, x = 1))+
    geom_quasirandom(alpha = 0.05)+
    geom_boxplot(outlier.shape = NA, fill = NA, color = "red")+
    scale_y_continuous(limits = c(0,4), breaks = seq(0,5, by = 0.5), position = "right")+
    theme_minimal()+
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank())+
    labs(x = "Whole Genome")

c + g + 
plot_layout(, nrow =1, widths = c(5,1))
