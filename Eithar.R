# set the working directory
setwd("E:/Either/new")

# read in read counts data and sample information
read_counts <- read.table("readcounts_featurecounts.tsv", sep = '\t', header = T, row.names = 1)
# read_counts data is read from a tab-separated file "readcounts_featurecounts.tsv"
# The first row is considered as the header and first column as row names

sample <- readxl::read_xls("sample.xls")
# sample information is read from an excel file "sample.xls" using the readxl library

# subset the sample information
sample <- sample[1:218,]
# The sample information is subsetted to only contain the first 218 rows
colnames(read_counts) <- sample$Run
# load the DESeq2 library
library(DESeq2)
# The DESeq2 library is loaded to perform differential gene expression analysis

# convert read counts to a matrix and remove missing values
read_counts <- as.matrix(read_counts)
read_counts <- na.omit(read_counts)
# The read counts data is converted to a matrix for further processing
# Missing values are removed using the "na.omit" function

# create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(read_counts, DataFrame(sample), ~ dead)
# A DESeqDataSet object is created from the read counts data and sample information
# The "~ dead" argument indicates that the "dead" variable in the sample information is used as a design formula for the analysis

# perform standard DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds)
# The standard DESeq2 analysis is performed on the DESeqDataSet object
# The results are stored in the "res" object

# summarize the results
summary(res)
# The results are summarized using the "summary" function

# perform moderated log2 fold changes analysis
resultsNames(dds)
resLFC <- lfcShrink(dds, coef=2, type="apeglm")
# Moderated log2 fold changes analysis is performed using the "lfcShrink" function
# The "coef" argument specifies the coefficient of interest (2) and the "type" argument specifies the type of shrinkage ("apeglm")
summary(resLFC)
# The results of the moderated log2 fold changes analysis are summarized

# store the differential expressed genes
DEGs <- as.data.frame(resLFC)
write.csv(DEGs, "Expression_values.csv")
# The differential expressed genes (DEGs) are stored as a data frame

# perform alternate analysis: likelihood ratio test
ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
resLRT <- results(ddsLRT)
# An alternate analysis is performed using the likelihood ratio test
# The "test" argument specifies the test to be used and the "reduced" argument specifies the reduced design formula (1)
# The results of the likelihood ratio test are stored in the "resLRT" object
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 10e-4,
                FCcutoff = 1.333,
                xlim = c(-5.5, 5.5),
                ylim = c(0, -log10(10e-12)),
                pointSize = 1.5,
                labSize = 2.5,
                title = 'DESeq2 results',
                subtitle = 'Differential expression',
                caption = 'FC cutoff, 1.333; p-value cutoff, 10e-4',
                legendPosition = "right",
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 0.9,
                drawConnectors = TRUE,
                hline = c(10e-8),
                widthConnectors = 0.5)
########################## PCA analysis ##########################################
vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:10]
df <- as.data.frame(colData(dds)[,c("dead","Smoking")])
DEGs_name <- rownames(DEGs)
pheatmap(assay(vsd)[DEGs_name[1:20],], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
plotPCA(vsd, intgroup=c("Smoking"))
pcaData <- plotPCA(vsd, intgroup=c("Smoking", "dead"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Smoking, shape=dead)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)[DEGs_name[1:20],]))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Smoking, vsd$dead, sep="-")
#colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
library(pheatmap)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
library(gplots)
heatmap.2(sampleDistMatrix, col=redgreen(250), scale="column", key=TRUE, keysize =0.9, symkey=FALSE, 
          density.info="none", trace="none", cexRow = 0.5, cexCol = 0.8, margins=c(5,7))
################## Functional Enrichment Analysis ################################
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl")
# host = "ensembl.org")
ttg <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id","external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = row.names(data),
  mart = mart)
library(topGO) #package for topGO gene ontology annotation
library(org.Hs.eg.db)
library("ggplot2")

data <- DEGs
data$ensembl_gene_id <- row.names(data)
data <- merge(x = data, y = ttg, by = "ensembl_gene_id", all = TRUE)
write.csv(data, "DEGs_with_geneName.csv")
#all <- read.csv("significant_de.csv", header = TRUE)
#View(data)
dim(data)
#dim(all)

up.idx <- which(data$pvalue < 0.05 & data$log2FoldChange > 1)
dn.idx <- which(data$pvalue < 0.05 & data$log2FoldChange < 1)

all.genes <- data$external_gene_name
up.genes <- data[up.idx,]$external_gene_name
dn.genes <- data[dn.idx,]$external_gene_name


ontology <- "CC" #"BP. CC, MF"
algorithm <- "weight01"
statistic <- "fisher"
outTitle <- paste0("topGO_GO-", ontology, "_ORA_", algorithm,"_", statistic)


upList <- factor(as.integer(all.genes %in% up.genes))
names(upList) <- all.genes
dnList <- factor(as.integer(all.genes %in% dn.genes))
names(dnList) <- all.genes



upGOdata <- new("topGOdata", ontology = ontology, allGenes = upList,geneSel = function(x)(x == 1),
                nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "SYMBOL")

dnGOdata <- new("topGOdata", ontology = ontology, allGenes = dnList,geneSel = function(x)(x == 1),
                nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "SYMBOL")

upRes <- runTest(upGOdata, algorithm = algorithm, statistic = statistic)
dnRes <- runTest(dnGOdata, algorithm = algorithm, statistic = statistic)


png(paste0(outTitle, "_up.png"), width = 8, height = 6, units = "in", res = 300)
enrichment_barplot(upGOdata, upRes, showTerms = 20, numChar = 50, orderBy = "Scores",
                   title = paste0("GO-", ontology,"up-regulated genes"))
invisible(dev.off())

png(paste0(outTitle, "_dn.png"), width = 8, height = 6, units = "in", res = 300)
enrichment_barplot(dnGOdata, dnRes, showTerms = 20, numChar = 50, orderBy = "Scores",
                   title = paste0("GO-", ontology," ORA of down-regulated genes"))
invisible(dev.off())


up.sigGenes <- sigGenes(upGOdata)
dn.sigGenes <- sigGenes(dnGOdata)
up.AnnoList <- lapply(up.tab$"GO.ID",
                      function(x) as.character(unlist(genesInTerm(object = upGOdata, whichGO = x))))
dn.AnnoList <- lapply(dn.tab$"GO.ID",
                      function(x) as.character(unlist(genesInTerm(object = dnGOdata, whichGO = x))))

up.tab <- GenTable(upGOdata, Pval = upRes, topNodes = 20)
dn.tab <- GenTable(dnGOdata, Pval = dnRes, topNodes = 20)

up.tab$Term <- sapply(up.tab$"GO.ID", function(go) Term(GO.db::GOTERM[[go]]))
dn.tab$Term <- sapply(dn.tab$"GO.ID", function(go) Term(GO.db::GOTERM[[go]]))

up.sigGenes <- sigGenes(upGOdata)
dn.sigGenes <- sigGenes(dnGOdata)

# Retrieve gene symbols for each GO from the test result
up.AnnoList <- lapply(up.tab$"GO.ID",
                      function(x) as.character(unlist(genesInTerm(object = upGOdata, whichGO = x))))
dn.AnnoList <- lapply(dn.tab$"GO.ID",
                      function(x) as.character(unlist(genesInTerm(object = dnGOdata, whichGO = x))))
up.SigList <- lapply(up.AnnoList, function(x) intersect(x, up.sigGenes))
dn.SigList <- lapply(dn.AnnoList, function(x) intersect(x, dn.sigGenes))
up.tab$Genes <- sapply(up.SigList, paste, collapse = ",")
dn.tab$Genes <- sapply(dn.SigList, paste, collapse = ",")
cbind(head(up.tab$Genes, 5))
cbind(head(dn.tab$Genes, 5))
write.table(up.tab, file = paste0(outTitle, "_up.txt"), sep = "\t", quote = F,
            row.names = F, col.names = T)
write.table(dn.tab, file = paste0(outTitle, "_dn.txt"), sep = "\t", quote = F,
            row.names = F, col.names = T)



#For KEGG Pathway Analysis
DisGeNET
Cancer_Cell_Line_Encyclopedia
ClinVar_2019
GO_Biological_Process_2021
GO_Cellular_Component_2021
GO_Molecular_Function_2021
Human_Phenotype_Ontology
KEGG_2019_Human
OMIM_Disease
Panther_2016
Reactome_2016
WikiPathways_2019_Human

library(enrichR)
dbs_pw <- "KEGG_2019_Human"
upEnriched_pw <- enrichr(genes = up.genes, databases = dbs_pw)
dnEnriched_pw <- enrichr(genes = dn.genes, databases = dbs_pw)

plotEnrich(upEnriched_pw[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value", title = "Up-regulated Genes - KEGG Pathway Analysis")
plotEnrich(dnEnriched_pw[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value", title = "Down-regulated Genes - KEGG Pathway Analysis")

write.csv(upEnriched_pw, "up_kegg NvR.csv")
write.csv(dnEnriched_pw, "dn_kegg NvR.csv")
KEGG_up <- upEnriched_pw$ClinVar_2019
allEnriched_pw <- enrichr(genes = all.genes, databases = dbs_pw)
View(allEnriched_pw$KEGG_2019_Human)
write.csv(allEnriched_pw$KEGG_2019_Human, "KEGG for NvR.csv")