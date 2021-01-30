library(SC3)
library(SingleCellExperiment)
library(cluster)
library(mclust)

# consensus approach
ARIs <- vector()
sils <- vector()
authors <- c('Biase', 'Goolam', 'Yan', 'Deng', 'Pollen', 'Bozec', 'Treutlein')
cells <- c(49, 124, 124, 259, 301, 1016, 201)
clusters <- c(3, 5, 8, 10, 11, 2, 4)
for (i in 1:length(authors)) {
	setwd(paste('../', sep = '', authors[i]))
	n <- cells[i]
	k <- clusters[i]
	ref <- scan('ref.txt')
	table <- read.table('norm.txt')
	dist <- as.dist(read.table('polaratio.txt'))
	# dist <- 1 - cor(table, method = 'spearman') # 'pearson'
	# dist <- dist(t(table)) # euclidean

	labels <- cutree(hclust(as.dist(dist)), k)
	ARIs <- append(ARIs, adjustedRandIndex(labels, ref))
	sils <- append(sils, mean(silhouette(labels, as.dist(dist))[, 3]))
	labels <- pam(as.dist(dist), k)[["clustering"]]
	ARIs <- append(ARIs, adjustedRandIndex(labels, ref))
	sils <- append(sils, mean(silhouette(labels, as.dist(dist))[, 3]))
	ARIs <- append(ARIs, 0)
	sils <- append(sils, -1)
	for (j in 1:5) {
		sc <- SingleCellExperiment(assays = list(counts = table, logcounts = as.matrix(table)))
		rowData(sc)$feature_symbol <- vector(mode = 'character', length = dim(sc@assays)[1])
		sc <- sc3_prepare(sc, gene_filter=F)
		sc@metadata[["sc3"]][["distances"]][["euclidean"]] <- as.matrix(dist)
		sc <- sc3_calc_consens(sc3_kmeans(sc3_calc_transfs(sc), k))
		labels <- reindex_clusters(sc@metadata$sc3$consensus[[as.character(k)]]$hc, k)
		ARI <- adjustedRandIndex(labels, ref)
		sil <- mean(silhouette(labels, dist)[, 3])
		if (sil > sils[length(sils)]) {
			sils[length(sils)] <- sil
			ARIs[length(ARIs)] <- ARI
		}
	}
	View(ARIs)
	View(sils)
}

# SC3
ARIs <- vector()
for (i in 1:length(authors)) {
	setwd(paste('../', sep = '', authors[i]))
	n <- cells[i]
	k <- clusters[i]
	ref <- scan('ref.txt')
	table <- read.table('norm.txt')

	for (j in 1:1) {
		sc <- SingleCellExperiment(assays = list(counts = table, logcounts = as.matrix(table)))
		rowData(sc)$feature_symbol <- vector(mode = 'character', length = dim(sc@assays)[1])

		sc <- sc3_calc_dists(sc3_prepare(sc, gene_filter = FALSE)) # default (euclidean, pearson, spearman)

		# individual distance metric
		# sc <- sc3_prepare(sc, gene_filter=F)
		# dist <- as.matrix(read.table('polaratio.txt'))
		# sc@metadata[["sc3"]][["distances"]][["euclidean"]] <- dist

		sc <- sc3_calc_consens(sc3_kmeans(sc3_calc_transfs(sc), k))
		labels <- reindex_clusters(sc@metadata$sc3$consensus[[as.character(k)]]$hc, k)
		ARIs <- append(ARIs, adjustedRandIndex(labels, ref))
	}
	View(ARIs)
}