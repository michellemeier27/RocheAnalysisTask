# Quality Control and Clean up of Data

# Set WD and load all packages ----
# set wd
setwd("//Users/michellemeier/RocheAnalysisTask/")
#libraries
library(GEOquery)
library(stringr)
library(PCAtools)
library(edgeR)
library(tidyverse)
library(reshape2)
library("circlize")
library("RColorBrewer")
library(ComplexHeatmap)



# Load data ----
CountMatrixSup <- read.delim('DataSet/GSE107593_raw_reads_BCHRNAseq.txt', check.names = FALSE)
C = colnames(CountMatrixSup)[10:57]
# load Meta Data with GEO reference
GDS <- getGEO(GEO = 'GSE107593',GSEMatrix = TRUE, getGPL = FALSE )


# Reformat Metadata ----
MetaDF = data.frame(SourceName = GDS$GSE107593_series_matrix.txt.gz$source_name_ch1,
                    Subject = GDS$GSE107593_series_matrix.txt.gz$`subject:ch1`,
                    InflammationStatus = str_replace(GDS$GSE107593_series_matrix.txt.gz$characteristics_ch1, 'status: ', ''),
                    Location = str_replace(GDS$GSE107593_series_matrix.txt.gz$characteristics_ch1.2, 'location: ', ''),
                    GSM = GDS$GSE107593_series_matrix.txt.gz$geo_accession)
#translate Names to match colnames in Countmatrix
MetaDF$SourceName = str_replace(MetaDF$SourceName, 'Colon_', '')
MetaDF$SourceName = str_replace(MetaDF$SourceName, ' ', '')
#spaces in C
C = str_replace(C, " ", '')
colnames(CountMatrixSup)[10:57] = C
#check if names are the same
all(sort(MetaDF$SourceName) == sort(C)) #yep, all is true
#reorder MetaDF
MetaDF = MetaDF[order(MetaDF$SourceName),]
rownames(MetaDF) = MetaDF$SourceName



# EdgeR ----
CountMatrix = CountMatrixSup[10:57];
rownames(CountMatrix) = CountMatrixSup$Row
#reorder to fit MetaData
CountMatrix = CountMatrix[order(colnames(CountMatrix))]
y <- DGEList(counts = CountMatrix, genes = CountMatrixSup[1:9], group = MetaDF$InflammationStatus)
keep <- filterByExpr(y) #filter out genes that are not expressed highly enough across all samples
y <- y[keep, , keep.lib.sizes=FALSE] # kicking them out + recalculating library size
y_norm <- calcNormFactors(y, method = 'TMM') #corrects for highly variable genes overshadowing all other effects
cpm <- cpm(y, log = TRUE) #CPM, without normalisation for highly variable genes
tmm <- cpm(y_norm, log = TRUE) #TMM, with normalisation for highly variable genes

# QC: boxplots ----
reshape_tmm <- melt(tmm)
colnames(reshape_tmm) = c('Gene', 'Sample', 'Count')
bp_tmm <- ggplot(reshape_tmm, aes(x=Sample, y=Count))+ geom_boxplot()+ggtitle('TMM normalised gene count distribution')
bp_tmm <- bp_tmm + theme(axis.text.x = element_text(angle = 90))
reshape_cpm <- melt(cpm)
colnames(reshape_cpm) = c('Gene', 'Sample', 'Count')
bp_cpm <- ggplot(reshape_cpm, aes(x=Sample, y=Count))+ geom_boxplot()+ggtitle('CPM normalised gene count distribution')
bp_cpm <- bp_cpm + theme(axis.text.x = element_text(angle = 90))
#save them
strPathSave = 'Results/QC'
ggsave('QC_TMMDistribution.png',plot = bp_tmm, path = strPathSave, device = 'png')
ggsave('QC_CPMDistribution.png',plot = bp_cpm, path = strPathSave, device = 'png')

# Heatmap for overview ----
#transfer normalised reads to zscore
zscore = scale(t(as.matrix(tmm)))
zscore = t(zscore)
type = MetaDF$Subject
ha = HeatmapAnnotation(
  df = data.frame(subject = type,
                  inflammation= MetaDF$InflammationStatus,
                  location = MetaDF$Location),
  annotation_height = unit(4, "mm")
)
t1 = Sys.time()
Heatmap(as.matrix(tmm), name = "z-score", top_annotation = ha,
        show_row_names = FALSE, show_column_names = FALSE)
#
 # Heatmap(CountMatrixSup$gene_type[keep], name = "Type", width = unit(5, "mm"))
t2 = Sys.time()
diff =t2-t1











