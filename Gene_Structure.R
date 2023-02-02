library(ggtranscript)
library(rtracklayer)
library(stringr)
library(ggplot2)
library(magrittr)
library(dplyr)
gff <- readGFF("G:/defensin_results/GCF_001704415.2_ARS1.2_genomic.gff", 
               tags=c("ID", "Parent", "Name"))
gff <- as.data.frame(gff)
gff$Name <- as.character(gff$Parent)
gff <- subset(gff, Parent %in% c("rna-XM_005674665.3"))
gff <- subset(gff, type == "exon")
gff %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = Name
  )) +
  geom_range(
    aes(fill = Name)
  ) +
  geom_intron(
    data = to_intron(gff, "Parent"),
    aes(strand = strand)
  )

                 