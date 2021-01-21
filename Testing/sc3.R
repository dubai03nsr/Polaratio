library(SC3)
library(SingleCellExperiment)

ARIs <- array(numeric())
table <- read.table('')
ref <- read.table('')
n <- 
k <- 
for (tt in 1:10) {
	sc <- SingleCellExperiment(assays = list(counts = table, logcounts = as.matrix(table)))
	rowData(sc)$feature_symbol <- vector(mode = 'character', length = dim(sc@assays)[1])
	sc <- sc3_calc_dists(sc3_prepare(sc, gene_filter=F))

	# align <- as.matrix(read.table('align.txt'))
	# polar <- as.matrix(read.table('polar.txt'))
	# sc@metadata[["sc3"]][["distances"]][["euclidean"]] <- align
	# sc@metadata[["sc3"]][["distances"]][["pearson"]] <- align
	# sc@metadata[["sc3"]][["distances"]][["spearman"]] <- polar

	sc <- sc3_calc_consens(sc3_kmeans(sc3_calc_transfs(sc), k))
	labels <- reindex_clusters(sc@metadata$sc3$consensus[[as.character(k)]]$hc, k)
	con <- matrix(0, nrow = k, ncol = k)
	for (i in 1:n) {
		con[labels[i], ref[i, 1] + 1] <- con[labels[i], ref[i, 1] + 1] + 1
	}
	a <- array(0, c(k))
	b <- array(0, c(k))
	for (i in 1:k) {
		for (j in 1:k) {
			a[i] <- a[i] + con[i, j]
			b[j] <- b[j] + con[i, j]
		}
	}
	s1 <- 0
	s2 <- 0
	s3 <- 0
	for (i in 1:k) {
		for (j in 1:k) {
			s1 <- s1 + con[i, j] * (con[i, j] - 1) / 2
			if (i == 1) s3 <- s3 + b[j] * (b[j] - 1) / 2
		}
		s2 <- s2 + a[i] * (a[i] - 1) / 2
	}
	nC2 <- n * (n - 1) / 2
	ARIs <- append(ARIs, (s1 - s2*s3/nC2) / (0.5 * (s2+s3) - s2*s3/nC2))
}
View(ARIs)