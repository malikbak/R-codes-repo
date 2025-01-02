library(rtracklayer)
gff <- readGFF("F:/Gulshan/TAIR10_GFF3_genes.gff")
gff <- as.data.frame(gff)
gff_tb <- read.table("F:/Gulshan/TAIR10_GFF3_genes.gff", comment.char = "#", sep = "\t")
vector <- unlist(strsplit(gff_tb$V9, ";"))
vector
library(stringr)
column <- unique(str_extract(vector, "[^=]+"))

res <- str_match(gff_tb$V9, sprintf("%s=\\s*(.*?)\\s*;", "ID"))
gff_tb["ID"] <- res[,2]
for (i in c(1:7)) {
  res <- str_match(gff_tb$V9, sprintf("%s=\\s*(.*?)\\s*;", column[i]))
  gff_tb[column[i]] <- res[,2]
}
