library(SC3)
library(SingleCellExperiment)
library(scater)

folder = ''
counts <- read.table(folder + '')
sc <- SingleCellExperiment(assays = list(counts = counts, logcounts = log2(counts + 1)))
rowData(sc)$feature_symbol <- vector(mode = 'character', length = dim(sc@assays)[1])
filt <- get_processed_dataset(sc3_prepare(sc))
norm <- logNormCounts(SingleCellExperiment(assays = list(counts = filt)))
write.table(norm@assays@data@listData[['logcounts']], file = folder + '', row.names=FALSE, col.names=FALSE)