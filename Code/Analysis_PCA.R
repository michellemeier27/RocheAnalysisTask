# PCA analysis

# Set WD and load all packages ----
# set wd
setwd("//Users/michellemeier/RocheAnalysisTask/")
#libraries
library(PCAtools)
library(edgeR)
library(tidyverse)
library(reshape2)
library("circlize")
library("RColorBrewer")
library(ComplexHeatmap)

# Run PCA ----
# keep the top 500 variable genes
top500 = 1- (500 / dim(CountMatrix)[1])
#use log transformed tmm (normal distribution assumed for PCA)
p <- pca(tmm, metadata = MetaDF, removeVar = top500)
screeplot(p, axisLabSize = 10, titleLabSize = 18) 
# scree plot shows that most variance is explained by PC1 (and PC2)

# Eigencorrelation ----
eigencorplot(p, metavars = colnames(MetaDF) , main = 'Eigencorrelation', cexMain = 1.5,
             cexLabX = 0.7, cexLabY = 0.7, cexLabColKey = 0.7, cexCorval = 0.7)
# strong correlation PC1 and Inflammation status => visualise in biplot

# Biplots ----
#colour only inflammation
biplot(p,
       lab = '',
       colby = 'InflammationStatus',
       hline = 0, vline = 0,
       legendPosition = 'right',
       title = 'PCA: Simple Overview',
       axisLabSize = 10)
strPathSave = 'Results/PCA'
ggsave('InflammationColbyBiplotPCA_log2TMM.png', path = strPathSave, device = 'png')
#color inflammation and label location
biplot(p,
       lab = p$metadata$Location,
       colby = 'InflammationStatus',
       hline = 0, vline = 0,
       legendPosition = 'right',
       title = 'PCA: Location + Inflammation',
       axisLabSize = 10)
strPathSave = 'Results/PCA'
ggsave('InflammationColbyLocationLabBiplotPCA_log2TMM.png', path = strPathSave, device = 'png')

#color location and label subject
biplot(p,
       lab = p$metadata$Subject,
       colby = 'Location',
       hline = 0, vline = 0,
       legendPosition = 'right',
       title = 'PCA: Location + Subject',
       axisLabSize = 10)
strPathSave = 'Results/PCA'
ggsave('LocationColbySubjectLabBiplotPCA_log2TMM.png', path = strPathSave, device = 'png')
# color location 
biplot(p,
       lab = '',
       colby = 'Location',
       hline = 0, vline = 0,
       legendPosition = 'right',
       title = 'PCA: Location',
       axisLabSize = 10)
strPathSave = 'Results/PCA'
ggsave('LocationColbybBiplotPCA_log2TMM.png', path = strPathSave, device = 'png')

#with loadings + inflammation colby
biplot(p,
       lab = '',
       colby = 'InflammationStatus',
       hline = 0, vline = 0,
       legendPosition = 'right',
       title = 'PCA: Loadings',
       axisLabSize = 10,
       showLoadings = TRUE,
       boxedLoadingsNames = FALSE)
strPathSave = 'Results/PCA'
ggsave('InflammationColbyLoadingsBiplotPCA_log2TMM.png', path = strPathSave, device = 'png')







p$loadings[1:5,1:5]
plotloadings(p, rangeRetain = 0.2) #make own plot
# Loadings ----
#access loading
Loadings = p$loadings
#PC1
SortPC1Loadings = Loadings[order(abs(Loadings[,1]), decreasing = TRUE),]
Top10PC1 = SortPC1Loadings[1:10,]
Top10PC1_ensemblnames = rownames(SortPC1Loadings)[1:10]
Top10PC1_normalnames = CountMatrixSup$gene_name[CountMatrixSup$Row %in% Top10PC1_ensemblnames]
Top10PC1_function = CountMatrixSup$gene_type[CountMatrixSup$Row %in% Top10PC1_ensemblnames]
melted_Loadings= melt(Loadings[,1:5])
highlighted = data.frame(variable = melted_Loadings$variable[1:10],
                         value = Top10PC1[,1],
                         name = Top10PC1_normalnames)
lp <- ggplot(melted_Loadings, aes(x=variable, y=value)) +
  geom_point() + geom_point(data = highlighted, aes(x=variable, y=value), color='red',
                            size=3) + geom_text_repel(data = highlighted, aes(x=variable, y=value, label=name),nudge_x = -0.6, max.overlaps = 50)
lp <- lp + ggtitle('Loadings PC1')
strPathSave = 'Results/PCA'
ggsave('LoadingsPC1PCA_log2TMM.png', path = strPathSave, device = 'png')

#PC2
SortPC2Loadings = Loadings[order(abs(Loadings[,2]), decreasing = TRUE),]
Top10PC2 = SortPC2Loadings[1:10,]
Top10PC2_ensemblnames = rownames(SortPC2Loadings)[1:10]
Top10PC2_normalnames = CountMatrixSup$gene_name[CountMatrixSup$Row %in% Top10PC2_ensemblnames]
Top10PC2_function = CountMatrixSup$gene_type[CountMatrixSup$Row %in% Top10PC2_ensemblnames]
melted_Loadings= melt(Loadings[,1:5])
highlighted = data.frame(variable = melted_Loadings$variable[129:138],
                         value = Top10PC2[,2],
                         name = Top10PC2_normalnames)
lp <- ggplot(melted_Loadings, aes(x=variable, y=value)) +
  geom_point() + geom_point(data = highlighted, aes(x=variable, y=value), color='red',
                            size=3) + geom_text_repel(data = highlighted, aes(x=variable, y=value, label=name),nudge_x = -0.6, max.overlaps = 50)
lp <- lp + ggtitle('Loadings PC2')
strPathSave = 'Results/PCA'
ggsave('LoadingsPC2PCA_log2TMM.png', path = strPathSave, device = 'png')

#summary df
SummaryLoadingsPC1 = cbind(Top10PC1[,1], Top10PC1_ensemblnames, Top10PC1_normalnames, Top10PC1_function)
colnames(SummaryLoadingsPC1) = c('Loading', 'ensembl', 'gene_names', 'gene_type')
write.csv(SummaryLoadingsPC1, 'Results/PCA/SummaryLoadingsPC1_log2TMM.csv')
SummaryLoadingsPC2 = cbind(Top10PC2[,2], Top10PC2_ensemblnames, Top10PC2_normalnames, Top10PC2_function)
colnames(SummaryLoadingsPC2) = c('Loading', 'ensembl', 'gene_names', 'gene_type')
write.csv(SummaryLoadingsPC2, 'Results/PCA/SummaryLoadingsPC2_log2TMM.csv')





