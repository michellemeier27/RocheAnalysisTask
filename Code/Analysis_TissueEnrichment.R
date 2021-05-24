# Tissue enrichment

# Set WD and load all packages ----
# set wd
setwd("//Users/michellemeier/RocheAnalysisTask/")
#libraries
library(TissueEnrich)

t <- decideTests(et, lfc = 1.5) #filter
rownames(t) = et$genes$gene_name
genes_1 <- rownames(t)[t == 1 | t == -1] #get genenames
gs <- GeneSet(geneIds = genes_1, organism = 'Homo Sapiens',  geneIdType= SymbolIdentifier())
output<-teEnrichment(inputGenes = gs)

seEnrichmentOutput<-output[[1]]
enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
enrichmentOutput$Tissue<-row.names(enrichmentOutput)
head(enrichmentOutput)

p <- ggplot(enrichmentOutput,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-log10(padj)')+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())
p <- p + ggtitle('Tissue Enrichment of significant DE (|log2FC| > 1.5)')
p
strPathSave = 'Results/TissueEnrichment'
ggsave('TissueEnrichment.png',plot = p, path = strPathSave, device = 'png')


