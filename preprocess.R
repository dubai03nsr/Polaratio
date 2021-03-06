if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
BiocManager::install("SC3")

library(SC3)
library(SingleCellExperiment)
library(scater)

folder = ''
counts <- read.table(paste(folder, sep='', 'counts.txt'))
sc <- SingleCellExperiment(assays = list(counts = counts, logcounts = log2(counts + 1)))
rowData(sc)$feature_symbol <- vector(mode = 'character', length = dim(sc@assays)[1])
filt <- get_processed_dataset(sc3_prepare(sc))
norm <- logNormCounts(SingleCellExperiment(assays = list(counts = filt)))
write.table(norm@assays@data@listData[['logcounts']], file = paste(folder, sep='', 'norm.txt'), row.names=FALSE, col.names=FALSE)