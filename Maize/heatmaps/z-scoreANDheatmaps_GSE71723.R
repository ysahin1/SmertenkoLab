library(NMF)
library(pheatmap)
zm_per <- read.delim("os_zm_sym_per_homologs.txt")
perox_GSE71723_de <- dn_GSE71723[which(dn_GSE71723 %in% zm_per$zm_locus)]
nperDEG <- length(perox_GSE71723_de)
nperDEG
sym_idx <- match(perox_GSE71723_de,zm_per$zm_locus)
sm_perox_GSE71723 <- as.vector(zm_per$symbols[sym_idx])
sm_perox_GSE71723
perox_vstc_GSE71723 <- vstc_dds_GSE71723[perox_GSE71723_de,]
dim(perox_vstc_GSE71723)
DEperox_zscore_GSE71723 <- (perox_vstc_GSE71723-rowMeans(perox_vstc_GSE71723))/(rowSds(as.matrix(perox_vstc_GSE71723)))
head(DEperox_zscore_GSE71723)
GSE71723_coln <- data.frame(sample = c("control","control","control","control","drought","drought","drought","drought"))
GSE71723_coln
row.names(GSE71723_coln) <- colnames(vstc_dds_GSE71723)
pheatmap(DEperox_zscore_GSE71723,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE71723_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE71723,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)

