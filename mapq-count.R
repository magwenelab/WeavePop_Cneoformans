library(tidyverse)
library(RColorBrewer)
library(scales)
library(svglite)
mapq<- read.delim("Downloads/mapq.tsv", header = FALSE, col.names = c("MAPQ", "Count"))

qual_colors <- colorRampPalette(brewer.pal(11, "RdYlBu"))(max(mapq$MAPQ)+1)

plot <- ggplot(mapq, aes(x=MAPQ, y=Count))+
  geom_col(fill = "hotpink4")+
  #scale_color_continuous(qual_colors)+
  #scale_fill_manual(qual_colors)+
  scale_y_continuous(name = "Count", labels = comma)+
  theme_light()
plot
ggsave("Downloads/mapq-count.svg", plot = plot, dpi = 200, units = "cm", height = 10, width = 10)
