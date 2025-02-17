# First look at data
#16.05.2021

# Set WD and load all packages ----
# set wd
setwd("//Users/michellemeier/RocheAnalysisTask/")
#libraries
library(GEOquery)
library(stringr)
library(PCAtools)
library(edgeR)

# Load data ----
CountMatrix <- read.delim('DataSet/GSE107593_raw_reads_BCHRNAseq.txt')
C = colnames(CountMatrix)[10:57]
# load Meta Data with GEO reference
GDS <- getGEO(GEO = 'GSE107593',GSEMatrix = TRUE, getGPL = FALSE )
gset <- getGEOSuppFiles('GSE107593')
fnames = rownames(gset)
b2 = read.delim(fnames[1], header = TRUE)

# Reformat Metadata ----
MetaDF = data.frame(SourceName = GDS$GSE107593_series_matrix.txt.gz$source_name_ch1,
                    Subject = GDS$GSE107593_series_matrix.txt.gz$`subject:ch1`,
                    InflammationStatus = str_replace(GDS$GSE107593_series_matrix.txt.gz$characteristics_ch1, 'status: ', ''),
                    Location = str_replace(GDS$GSE107593_series_matrix.txt.gz$characteristics_ch1.2, 'location: ', ''),
                    GSM = GDS$GSE107593_series_matrix.txt.gz$geo_accession)
#translate Names to match colnames in Countmatrix
translatedNames = str_replace(MetaDF$SourceName, 'Colon_', 'X')
translatedNames = str_replace(translatedNames, "[+-]", '.')
translatedNames = str_replace(translatedNames, " ", '')
#spaces in C
C = str_replace(translatedNames, " ", '')
#check if names are the same
all(sort(translatedNames) == sort(C)) #yep, all is true
#add to MetaDF
MetaDF$TranslatedNames = translatedNames

# EdgeR
CM = CountMatrix[10:57];
rownames(CM) = CountMatrix$Row
y <- DGEList(counts = CM)
keep <- filterByExpr(y)
y <- calcNormFactors(y, method = 'TMM')
cpm <- cpm(y, log = TRUE)

#PCA
p <- pca(cpm, removeVar = 0.9917353)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p)
p$loadings[1:5,1:5]
plotloadings(p, rangeRetain = 0.2) #make own plot












