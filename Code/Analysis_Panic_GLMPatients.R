# DE analysis v2.1
# difference to version 2.0: this is with linear model 

# Slrt WD and load all packages ----
# slrt wd 
slrtwd("//Users/michellemeier/RocheAnalysisTask/")
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


#glm model ----
patient  <- as.factor(MlrtaDF$Subject)
location <- as.factor(MlrtaDF$Location)
iF <- MlrtaDF$InflammationStatus == 'Non-Inflamed'
iT <- MlrtaDF$InflammationStatus == 'Inflamed'
inflammation <- rep('T', times = dim(MlrtaDF)[1])
inflammation[iF] = 'F' 
inflammation <- as.factor(inflammation)
design_patients <- model.matrix(~patient+inflammation, data = y$samples)
rownames(design_patients) <- colnames(y)

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



# EdgeR: with design glm ----
#biological variation
y <- estimateDisp(y, design_patients, robust = TRUE)
BV <- sqrt(y$common.dispersion)
plotBCV(y) # read up on this, i don't really know what this means
#DE analysis
fit <- glmFit(y, design_patients)
lrt <- glmLRT(fit)


TopHits20Pval = topTags(lrt, n = 20) # internlrt says filter with pvals
write.csv(TopHits20Pval$table, 'Results/DE_GLMPatients/SummaryTop20Hits_GLMPatiens.csv')
#summary of up and down regulated genes
SummaryPval = summary(decideTests(lrt))
SummaryFCandPval = summary(decideTests(lrt, lfc = 1.5)) #apparently not recommended to also have fc cutoff
#make pie charts
pie(SummaryPval, rownames(SummaryPval), main = 'Summary DE analysis: filter FDR < 0.05')
pie(SummaryFCandPval, rownames(SummaryFCandPval), main = 'Summary DE analysis: filter FDR < 0.05 & log2FC > 1.5')


#volcano plot
df_plotting = data.frame(log2FC = lrt$table$logFC, FDR = lrt$table$PValue)
g <- ggplot(data = df_plotting, aes(x = log2FC, y = -log10(FDR))) +geom_point() 
lrt_low = data.frame(FC = lrt$table$logFC[lrt$table$logFC < -1.5 & lrt$table$PValue < 0.05],
                     FDR = lrt$table$PValue[lrt$table$logFC < -1.5 & lrt$table$PValue < 0.05])
lrt_high = data.frame(FC = lrt$table$logFC[lrt$table$logFC > 1.5 & lrt$table$PValue < 0.05],
                      FDR = lrt$table$PValue[lrt$table$logFC > 1.5 & lrt$table$PValue < 0.05])
g <- g + geom_point(data = lrt_low, aes(x = FC, y = -log10(FDR)), color='red')
g <- g + geom_point(data = lrt_high, aes(x = FC, y = -log10(FDR)), color='blue')
g <- g +ggtitle('Volcano Plot: FDR < 0.05, |log2FC| > 1.5')
strPathSave = 'Results/DE_GLMPatients'
ggsave('DE_VolcanoPlotFDRlFC.png',plot = g, path = strPathSave, device = 'png')

#volcano plot
df_plotting = data.frame(log2FC = lrt$table$logFC, FDR = lrt$table$PValue)
g <- ggplot(data = df_plotting, aes(x = log2FC, y = -log10(FDR))) +geom_point() 
lrt_low = data.frame(FC = lrt$table$logFC[lrt$table$logFC < 0 & lrt$table$PValue < 0.05],
                     FDR = lrt$table$PValue[lrt$table$logFC < 0& lrt$table$PValue < 0.05])
lrt_high = data.frame(FC = lrt$table$logFC[lrt$table$logFC > 0& lrt$table$PValue < 0.05],
                      FDR = lrt$table$PValue[lrt$table$logFC >0 & lrt$table$PValue < 0.05])
g <- g + geom_point(data = lrt_low, aes(x = FC, y = -log10(FDR)), color='red')
g <- g + geom_point(data = lrt_high, aes(x = FC, y = -log10(FDR)), color='blue')
g <- g +ggtitle('Volcano Plot: FDR < 0.05')
strPathSave = 'Results/DE_GLMPatients'
ggsave('DE_VolcanoPlotFDR.png',plot = g, path = strPathSave, device = 'png')


# Reactome stuff ----
genelist = -log10(lrt$table$PValue)*lrt$table$logFC
names(genelist) = lrt$genes$Entrez
genelist <- sort(genelist, decreasing = TRUE)
ge = gsePathway(genelist,  eps = 0)
PathwayList = ge@result$Description
#plot-options
cp <- cnetplot(ge, node_label="category") 
cp <- cp + ggtitle('Pathway Enrichment (GSEA): Pre-ranked -log10(pval) * log2FC')
ggsave('Reactome_GSEA_cnlrtplot_pvalFC.png',plot = cp, path = strPathSave, device = 'png')
dp <- dotplot(ge, showCategory=15)
dp <- dp + ggtitle('Pathway Enrichment (GSEA): Pre-ranked -log10(pval) * log2FC')
ggsave('Reactome_GSEA_dotplot_pvalFC.png',plot = dp, path = strPathSave, device = 'png')




