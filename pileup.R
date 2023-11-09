library(tidyverse)
library(RColorBrewer)
library(scales)
library(svglite)
mapq_full<- read.delim("Downloads/snps.pileup", header = FALSE, col.names = c("Chromosome", "Position", "Reference_Base", "Read_Count", "Read_Results", "Quality", "MAPQ" ))

mapq <- mapq_full %>%
  filter(Read_Count != 0)%>%
  filter(Read_Count < 10)%>%
  select(Chromosome, Position, Read_Count, MAPQ)

mat <- as.data.frame(str_split_fixed(mapq$MAPQ, ",", n = max(mapq$Read_Count)))

mapq <- mapq %>%
  select(-MAPQ)

mapq <- bind_cols(mapq, mat)

mat <- mat %>%
  rowwise()%>%
  mutate(Mean = mean(select(.,everything())))
