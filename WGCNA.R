# Uncomment and modify the following to install any missing packages
# install.packages(c("tidyverse", "magrittr", "WGCNA))
library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)      # provides the %>% operator
library(WGCNA)
library(DESeq2)
vsd <- varianceStabilizingTransformation(dds)
library(genefilter)
wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]
write.csv(expr_normalized, "Normalized expression.csv")
expr_normalized[1:5,1:10]
expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

association <- subset(expr_normalized_df, Gene_id %in% data$ensembl_gene_id)
sample$name <- sample$Run
condi <- sample[,c(31,34,36)]
association <- merge(x = association, y = condi, by = "name", all.x = TRUE)
data$Gene_id <- data$ensembl_gene_id
association <- merge(x = association, y = data, by = "Gene_id", all.x = TRUE)
association$dead[association$dead == 0] <- "dead"
association$dead[association$dead == 1] <- "Live"
write.csv(association, "Association of gene with stages.csv")
up_associa <- subset(association, log2FoldChange > 1.5)
up_associa <- subset(up_associa, pvalue < 0.05)

col = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
                 "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                 "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
                 "#8A7C64", "#599861")
library(ggplot2)
ggplot(up_associa, aes(fill = dead,y=value, x=external_gene_name)) + 
  geom_bar(position="dodge", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=sample(col,2))
dn_associa <- subset(association, log2FoldChange < 1.5)
dn_associa <- subset(dn_associa, pvalue < 0.05)
ggplot(dn_associa, aes(fill = dead,y=value, x=external_gene_name)) + 
  geom_bar(position="dodge", stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=sample(col,2))
input_mat = t(expr_normalized)

input_mat[1:5,1:10]
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")
picked_power = 9
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor 
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )
#Relating modules to characteristics and identifying important genes
#Defining the number of genes and samples
nGenes = ncol(input_mat)
nSamples = nrow(input_mat)

#Recalculating MEs with label colors
moduleColors = mergedColors
MEs0 = moduleEigengenes(input_mat, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
datTraits = sample[,c(31,34)]
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
#sizeGrWindow(8,4)

#Displaying correlations and its p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

#Displaying the correlation values in a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
Peso10dias = as.data.frame(datTraits$Smoking)
names(Peso10dias) = "Smoking"
geneTraitSignificance = as.data.frame(cor(input_mat, Peso10dias, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
geneInfo0 = data.frame(ensembl_gene_id = colnames(input_mat),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
geneInfo0 <- merge(x = geneInfo0, y = data, by = "ensembl_gene_id", all.x = TRUE)
write.csv(geneInfo0, "Correlation_smoking.csv")
######################### DEAD ######################################
Peso10dias = as.data.frame(datTraits$dead)
names(Peso10dias) = "Dead"
geneTraitSignificance = as.data.frame(cor(input_mat, Peso10dias, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
geneInfo0 = data.frame(ensembl_gene_id = colnames(input_mat),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
geneInfo0 <- merge(x = geneInfo0, y = data, by = "ensembl_gene_id", all.x = TRUE)
write.csv(geneInfo0, "Correlation_dead.csv")
