# Pathway Analyis (Reactome)

# Set WD and load all packages ----
# set wd
setwd("//Users/michellemeier/RocheAnalysisTask/")
#libraries
library(edgeR)
library(tidyverse)
library('statmod')
library("AnnotationDbi")
library(org.Hs.eg.db)
library(ChIPpeakAnno)
library(stringr)
library(ReactomePA)
library(signatureSearch)
library('ggnewscale')
library(clusterProfiler)
library('ggupset')
strPathSave = 'Results/Reactome'

# pathway enrichment (reactome) for DE genes----
# all sig. genes abs(fc) > 1.5
t <- decideTests(et, lfc = 1.5) #filter
rownames(t) = et$genes$Entrez
genes <- rownames(t)[t == 1 | t == -1] #get genenames
en = enrichPathway(gene= genes, readable = T)
dp <- dotplot(en, showCategory=15) 
dp <- dp + ggtitle('ORA: |log2FC| > 1.5')
ggsave('Reactome_ORA_dotplot.png',plot = dp, path = strPathSave, device = 'png')
cp <- cnetplot(en, node_label="category") 
cp <- cp + ggtitle('ORA: |log2FC| > 1.5')
ggsave('Reactome_ORA_cnetplot.png',plot = cp, path = strPathSave, device = 'png')


# Reactome GSEA enrichment ----
#make gene list pvalue and fc interaction score
genelist = -log10(et$table$PValue)*et$table$logFC
names(genelist) = et$genes$Entrez
genelist <- sort(genelist, decreasing = TRUE)
ge = gsePathway(genelist,  eps = 0)
PathwayList = ge@result$Description
write.csv(PathwayList, 'Results/DE/SignificantPathwaysGSEA_pvalFC.csv')
#plot-options
cp <- cnetplot(ge, node_label="category") 
cp <- cp + ggtitle('Pathway Enrichment (GSEA): Pre-ranked -log10(pval) * log2FC')
ggsave('Reactome_GSEA_cnetplot_pvalFC.png',plot = cp, path = strPathSave, device = 'png')
dp <- dotplot(ge, showCategory=15)
dp <- dp + ggtitle('Pathway Enrichment (GSEA): Pre-ranked -log10(pval) * log2FC')
ggsave('Reactome_GSEA_dotplot_pvalFC.png',plot = dp, path = strPathSave, device = 'png')
up <- upsetplot(ge)
up <- up + ggtitle('Pathway Enrichment (GSEA): Pre-ranked -log10(pval) * log2FC')
#emapplot(pairwise_termsim(ge)) # don't really like this plot
vp <- viewPathway('Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell', readable = T, foldChange = genelist)
vp <- vp + ggtitle('Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell')
ggsave('Reactome_GSEA_vpImmunoreg_pvalFC.png',plot = vp, path = strPathSave, device = 'png')

#make gene list pvalue
genelist = -log10(et$table$PValue)
names(genelist) = et$genes$Entrez
genelist <- sort(genelist, decreasing = TRUE)
ge = gsePathway(genelist,  eps = 0)
PathwayList = ge@result$Description
write.csv(PathwayList, 'Results/DE/SignificantPathwaysGSEA_pval.csv')
#plot-options
cp <- cnetplot(ge, node_label="category") 
cp <- cp + ggtitle('Pathway Enrichment (GSEA): Pre-ranked -log10(pval)')
ggsave('Reactome_GSEA_cnetplot_pval.png',plot = cp, path = strPathSave, device = 'png')
dp <- dotplot(ge, showCategory=15)
dp <- dp + ggtitle('Pathway Enrichment (GSEA): Pre-ranked -log10(pval)')
ggsave('Reactome_GSEA_dotplot_pval.png',plot = dp, path = strPathSave, device = 'png')
up <- upsetplot(ge)
up <- up + ggtitle('Pathway Enrichment (GSEA): Pre-ranked -log10(pval)')
ggsave('Reactome_GSEA_upsetplot_pval.png',plot = up, path = strPathSave, device = 'png')

#emapplot(pairwise_termsim(ge)) # don't really like this plot
vp <- viewPathway('Integrin cell surface interactions', readable = T, foldChange = genelist)
vp <- vp + ggtitle('Integrin cell surface interactions')
ggsave('Reactome_GSEA_vpIntegrin_pval.png',plot = vp, path = strPathSave, device = 'png')







