library(SingleCellExperiment)
library(SC3)
library(cluster)

n <-  # number of cells
k <-  # number of clusters
table <- read.table('.txt') # expression matrix
dist <- as.dist(read.table('PD.txt')) # distance matrix
maxSil <- -1
select <- vector()

HC <- cutree(hclust(dist), k)
HCsil <- mean(silhouette(HC, dist)[, 3])
if (HCsil > maxSil) {
	maxSil <- HCsil
	select <- HC
}

KM <- pam(as.dist(dist), k)[["clustering"]]
KMsil <- mean(silhouette(KM, dist)[, 3])
if (KMsil > maxSil) {
	maxSil <- KMsil
	select <- KM
}

for (run in 1:1) {
	sc <- SingleCellExperiment(assays = list(counts = table, logcounts = as.matrix(table)))
	rowData(sc)$feature_symbol <- vector(mode = 'character', length = dim(sc@assays)[1])
	sc <- sc3_prepare(sc, gene_filter = FALSE)
	sc@metadata[["sc3"]][["distances"]][["euclidean"]] <- as.matrix(dist)
	sc <- sc3_calc_consens(sc3_kmeans(sc3_calc_transfs(sc), k))
	SC3 <- reindex_clusters(sc@metadata$sc3$consensus[[as.character(k)]]$hc, k)
	SC3sil <- mean(silhouette(SC3, dist)[, 3])
	if (SC3sil > maxSil) {
		maxSil <- SC3sil
		select <- SC3
	}
}

write.table(select, 'clustering.txt', row.names = FALSE, col.names = FALSE) # output file