# PCA analysis v2.0

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
library(PCAtools)
library(edgeR)
library(tidyverse)
library('statmod')
library("AnnotationDbi")
library(org.Hs.eg.db)
library(ChIPpeakAnno)
library(stringr)
library(ReactomePA)
library('biomaRt')


# adapting to ensembl gene id issue ----
#new approach
idfound <- y$genes$gene_name %in% mappedRkeys(org.Hs.egALIAS2EG)
y <- y[idfound,]
egSYMBOL <- toTable(org.Hs.egSYMBOL)
egEntrez <- toTable(org.Hs.egALIAS2EG)
m <- match(y$genes$gene_name, egEntrez$alias_symbol)
y$genes$Entrez <- egEntrez$gene_id[m]

o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$Entrez)
y <- y[!d,]
nrow(y)
#recomputing lib size
y$samples$lib.size <- colSums(y$counts)
#do cpm and tmm calculations again
y_norm <- calcNormFactors(y, method = 'TMM') #corrects for highly variable genes overshadowing all other effects
cpm <- cpm(y, log = TRUE) #CPM, without normalisation for highly variable genes
tmm <- cpm(y_norm, log = TRUE) #TMM, with normalisation for highly variable genes

# i don't like ensembl gene names
gnames <- rownames(tmm)
gsymbols <- CountMatrixSup$gene_name[match(gnames, CountMatrixSup$Row)]

tmm_symbols = tmm
rownames(tmm_symbols) = gsymbols
# Run PCA ----
set.seed(2021)
# keep the top 500 variable genes
top500 = 1- (500 / dim(CountMatrix)[1])
#use log transformed tmm (normal distribution assumed for PCA)
p <- pca(tmm, metadata = MetaDF, removeVar = top500)
p_new <- pca(tmm_symbols, metadata = MetaDF, removeVar = top500)
screeplot(p, axisLabSize = 10, titleLabSize = 18) 
# scree plot shows that most variance is explained by PC1 (and PC2)
strPathSave = 'Results/PCA_v2.0'
ggsave('ScreePlot_log2TMM.png', path = strPathSave, device = 'png')


# Eigencorrelation ----
eigencorplot(p, metavars = colnames(MetaDF) , main = 'Eigencorrelation', cexMain = 1.5,
             cexLabX = 0.7, cexLabY = 0.7, cexLabColKey = 0.7, cexCorval = 0.7)
# strong correlation PC1 and Inflammation status => visualise in biplot
strPathSave = 'Results/PCA_v2.0'
ggsave('Eigencorrelation_log2TMM.png', path = strPathSave, device = 'png')

# Biplots ----
#colour only inflammation
biplot(p,
       lab = '',
       colby = 'InflammationStatus',
       hline = 0, vline = 0,
       legendPosition = 'right',
       title = 'PCA: Simple Overview',
       axisLabSize = 10)
ggsave('InflammationColbyBiplotPCA_log2TMM.png', path = strPathSave, device = 'png')

#color inflammation and label location
biplot(p,
       lab = '',
       colby = 'Subject',
       hline = 0, vline = 0,
       legendPosition = 'right',
       title = 'PCA: Subject',
       axisLabSize = 10)
ggsave('SubjectColbyBiplotPCA_log2TMM.png', path = strPathSave, device = 'png')

#color location and label subject
biplot(p,
       lab = '',
       colby = 'Location',
       hline = 0, vline = 0,
       legendPosition = 'right',
       title = 'PCA: Location',
       axisLabSize = 10)
ggsave('LocationColbyBiplotPCA_log2TMM.png', path = strPathSave, device = 'png')

#with loadings + inflammation colby
biplot(p_new,
       lab = '',
       colby = 'InflammationStatus',
       hline = 0, vline = 0,
       legendPosition = 'right',
       title = 'PCA: Loadings',
       axisLabSize = 10,
       showLoadings = TRUE,
       boxedLoadingsNames = FALSE)
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
Top10PC1_function = CountMatrixSup$chr[CountMatrixSup$Row %in% Top10PC1_ensemblnames]
melted_Loadings= melt(Loadings[,1:5])
highlighted = data.frame(variable = melted_Loadings$variable[1:10],
                         value = Top10PC1[,1],
                         name = Top10PC1_normalnames)
lp <- ggplot(melted_Loadings, aes(x=variable, y=value)) +
  geom_point() + geom_point(data = highlighted, aes(x=variable, y=value), color='red',
                            size=3) + geom_text_repel(data = highlighted, aes(x=variable, y=value, label=name),nudge_x = -0.6, max.overlaps = 50)
lp <- lp + ggtitle('Loadings PC1')
ggsave('LoadingsPC1PCA_log2TMM.png', path = strPathSave, device = 'png')

#PC2
SortPC2Loadings = Loadings[order(abs(Loadings[,2]), decreasing = TRUE),]
Top10PC2 = SortPC2Loadings[1:10,]
Top10PC2_ensemblnames = rownames(SortPC2Loadings)[1:10]
Top10PC2_normalnames = CountMatrixSup$gene_name[match(Top10PC2_ensemblnames, CountMatrixSup$Row)]
melted_Loadings= melt(Loadings[,1:5])
highlighted = data.frame(variable = melted_Loadings$variable[129:138],
                         value = Top10PC2[,2],
                         name = Top10PC2_normalnames)
lp <- ggplot(melted_Loadings, aes(x=variable, y=value)) +
  geom_point() + geom_point(data = highlighted, aes(x=variable, y=value), color='red',
                            size=3) + geom_text_repel(data = highlighted, aes(x=variable, y=value, label=name),nudge_x = -0.6, max.overlaps = 50)
lp <- lp + ggtitle('Loadings PC2')
ggsave('LoadingsPC2PCA_log2TMM.png', path = strPathSave, device = 'png')

#summary df
SummaryLoadingsPC1 = cbind(Top10PC1[,1], Top10PC1_ensemblnames, Top10PC1_normalnames, Top10PC1_function)
colnames(SummaryLoadingsPC1) = c('Loading', 'ensembl', 'gene_names', 'gene_type')
write.csv(SummaryLoadingsPC1, 'Results/PCA_v2.0/SummaryLoadingsPC1_log2TMM.csv')
SummaryLoadingsPC2 = cbind(Top10PC2[,2], Top10PC2_ensemblnames, Top10PC2_normalnames, Top10PC2_function)
colnames(SummaryLoadingsPC2) = c('Loading', 'ensembl', 'gene_names', 'gene_type')
write.csv(SummaryLoadingsPC2, 'Results/PCA_v2.0/SummaryLoadingsPC2_log2TMM.csv')

# rest
pairsplot(p)



