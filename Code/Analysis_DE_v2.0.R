# DE analysis v2.0
# difference to first version: made the genes entrez compatible for GSEA after
# if one gene symbol id matched to multiple entrez, the first one was used 
# if multiple gene symbol id matched to the same entrez, the most highly expressed one was used 

# Set WD and load all packages ----
# set wd 
setwd("//Users/michellemeier/RocheAnalysisTask/")
#libraries
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




# EdgeR: without design glm ----
#biological variation
y <- estimateDisp(y, robust = TRUE)
BV <- sqrt(y$common.dispersion)
plotBCV(y) # read up on this, i don't really know what this means
#DE analysis
et <- exactTest(y, pair = c('Non-Inflamed', 'Inflamed')) #change order to see changes in inflamed
TopHits20Pval = topTags(et, n = 20) # internet says filter with pvals
write.csv(TopHits20Pval$table, 'Results/DE/SummaryTop20Hits.csv')
#summary of up and down regulated genes
SummaryPval = summary(decideTests(et))
SummaryFCandPval = summary(decideTests(et, lfc = 1.5)) #apparently not recommended to also have fc cutoff
#make pie charts
pie(SummaryPval, rownames(SummaryPval), main = 'Summary DE analysis: filter FDR < 0.05')
pie(SummaryFCandPval, rownames(SummaryFCandPval), main = 'Summary DE analysis: filter FDR < 0.05 & log2FC > 1.5')

#check correlation FC and p vals => if correlation then makes sense to have fc cutoff
ggplot(et$table, aes(x =logFC^2, y = PValue)) + geom_point()
cor(et$table$PValue, et$table$logFC^2) #no correlation

#volcano plot
df_plotting = data.frame(log2FC = et$table$logFC, FDR = et$table$PValue)
g <- ggplot(data = df_plotting, aes(x = log2FC, y = -log10(FDR))) +geom_point() 
et_low = data.frame(FC = et$table$logFC[et$table$logFC < -1.5 & et$table$PValue < 0.05],
                    FDR = et$table$PValue[et$table$logFC < -1.5 & et$table$PValue < 0.05])
et_high = data.frame(FC = et$table$logFC[et$table$logFC > 1.5 & et$table$PValue < 0.05],
                     FDR = et$table$PValue[et$table$logFC > 1.5 & et$table$PValue < 0.05])
g <- g + geom_point(data = et_low, aes(x = FC, y = -log10(FDR)), color='red')
g <- g + geom_point(data = et_high, aes(x = FC, y = -log10(FDR)), color='blue')
g <- g +ggtitle('Volcano Plot: FDR < 0.05, |log2FC| > 1.5')
strPathSave = 'Results/DE'
ggsave('DE_VolcanoPlotFDRlFC.png',plot = g, path = strPathSave, device = 'png')

#volcano plot
df_plotting = data.frame(log2FC = et$table$logFC, FDR = et$table$PValue)
g <- ggplot(data = df_plotting, aes(x = log2FC, y = -log10(FDR))) +geom_point() 
et_low = data.frame(FC = et$table$logFC[et$table$logFC < 0 & et$table$PValue < 0.05],
                    FDR = et$table$PValue[et$table$logFC < 0& et$table$PValue < 0.05])
et_high = data.frame(FC = et$table$logFC[et$table$logFC > 0& et$table$PValue < 0.05],
                     FDR = et$table$PValue[et$table$logFC >0 & et$table$PValue < 0.05])
g <- g + geom_point(data = et_low, aes(x = FC, y = -log10(FDR)), color='red')
g <- g + geom_point(data = et_high, aes(x = FC, y = -log10(FDR)), color='blue')
g <- g +ggtitle('Volcano Plot: FDR < 0.05')
strPathSave = 'Results/DE'
ggsave('DE_VolcanoPlotFDR.png',plot = g, path = strPathSave, device = 'png')




# Is the expression profile the same for top20hits across patients? ----
#only FDR hits
ix_ordered = order(et$table$PValue)
Top20AcrossPatients = tmm[ix_ordered[1:20],]
new_rowlabels = CountMatrixSup$gene_name[match(rownames(Top20AcrossPatients) , CountMatrixSup$Row)]
chromosome_location = et$genes$chr[ix_ordered[1:20]]
# make heatmap for visualisation
mat = scale(t(as.matrix(Top20AcrossPatients)))
mat = t(mat)
type = MetaDF$Subject
ha = HeatmapAnnotation(
  df = data.frame(subject = type,
                  inflammation= MetaDF$InflammationStatus),
  annotation_height = unit(4, "mm"))
Heatmap(mat, name = "z-score", km = 4, top_annotation = ha, 
        row_labels = new_rowlabels, column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 8))


#DFR+FC
et_table_subset = subset(et$table, abs(et$table$logFC) > 1.5)
ix_ordered = rownames(et_table_subset)[order(et_table_subset$PValue)]
Top20AcrossPatients = tmm[match(ix_ordered[1:20], rownames(tmm)),]
new_rowlabels = CountMatrixSup$gene_name[match(rownames(Top20AcrossPatients) , CountMatrixSup$Row)]
# make heatmap for visualisation
mat = scale(t(as.matrix(Top20AcrossPatients)))
mat = t(mat)
type = MetaDF$Subject
ha = HeatmapAnnotation(
  df = data.frame(subject = type,
                  inflammation= MetaDF$InflammationStatus),
  annotation_height = unit(4, "mm"))
Heatmap(mat, name = "zscore", km = 4, top_annotation = ha, 
        row_labels = new_rowlabels, column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 8))
