# DE analysis

# Set WD and load all packages ----
# set wd
setwd("//Users/michellemeier/RocheAnalysisTask/")
#libraries
library(PCAtools)
library(edgeR)
library(tidyverse)
library('statmod')


# EdgeR: without design glm ----
y <- estimateDisp(y, robust = TRUE)
BV <- sqrt(y$common.dispersion)
plotBCV(y) # read up on this, i don't really know what this means
et <- exactTest(y, pair = c('Non-Inflamed', 'Inflamed')) #change order to see changes in inflamed
TopHits20Pval = topTags(et, n = 20) # internet says filter with pvals
write.csv(TopHits20Pval$table, 'Results/DE/SummaryTop20Hits.csv')
#summary
SummaryPval = summary(decideTests(et))
SummaryFCandPval = summary(decideTests(et, lfc = 1.5)) #apparently not recommended to also have fc cutoff
#make pie charts
pie(SummaryPval, rownames(SummaryPval), main = 'Summary DE analysis: filter FDR < 0.05')
pie(SummaryFCandPval, rownames(SummaryFCandPval), main = 'Summary DE analysis: filter FDR < 0.05 & log2FC > 1.5')
#check correlation FC and p vals
ggplot(df_plotting, aes(x =log2FC^2, y = pval)) + geom_point()
cor(df_plotting$pval, df_plotting$log2FC^2) #no correlation
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


# Is the expression profile the same for top20hits across patients? ----
#only FDR hits
ix_ordered = order(et$table$PValue)
Top20AcrossPatients = tmm[ix_ordered[1:20],]
# make heatmap for visualisation
mat = as.matrix(Top20AcrossPatients)
type = MetaDF$Subject
ha = HeatmapAnnotation(
  df = data.frame(subject = type,
                  inflammation= MetaDF$InflammationStatus),
  annotation_height = unit(4, "mm"))
Heatmap(mat, name = "expression", km = 4, top_annotation = ha,
        show_row_names = FALSE, column_names_gp = gpar(fontsize = 8))

#DFR+FC
et_table_subset = subset(et$table, abs(et$table$logFC) > 1.5)
ix_ordered = order(et_table_subset$PValue)
Top20AcrossPatients = tmm[ix_ordered[1:20],]
# make heatmap for visualisation
mat = as.matrix(Top20AcrossPatients)
type = MetaDF$Subject
ha = HeatmapAnnotation(
  df = data.frame(subject = type,
                  inflammation= MetaDF$InflammationStatus),
  annotation_height = unit(4, "mm"))
Heatmap(mat, name = "expression", km = 4, top_annotation = ha,
        show_row_names = FALSE, column_names_gp = gpar(fontsize = 8))



# Edge R with design glm ----
# although i don't think we necessarily need it...
#making design matrix
Subject <- as.factor(MetaDF$Subject)
Location <- as.factor(MetaDF$Location)
InflammationStatus <- as.factor(MetaDF$InflammationStatus)
design <- model.matrix(~ Subject + Location + InflammationStatus)
rownames(design) = colnames(y)
# run
y <- estimateDisp(y, design = design, robust = TRUE)
BV <- sqrt(y$common.dispersion)
plotBCV(y) # read up on this, i don't really know what this means
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)
df_plotting = data.frame(log2FC = lrt$table$logFC, pval = lrt$table$PValue)
ggplot(data = df_plotting, aes(x = log2FC, y = -log10(pval))) +geom_point()



