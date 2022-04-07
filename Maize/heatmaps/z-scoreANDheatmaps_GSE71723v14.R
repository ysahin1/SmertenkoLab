library(NMF)
library(pheatmap)
zm_per <- read.delim("os_zm_sym_per_homologs.txt")
perox_GSE71723v14_de <- dn_GSE71723v14[which(dn_GSE71723v14 %in% zm_per$zm_locus)]
nperDEG <- length(perox_GSE71723v14_de)
nperDEG
sym_idx <- match(perox_GSE71723v14_de,zm_per$zm_locus)
sm_perox_GSE71723v14 <- as.vector(zm_per$symbols[sym_idx])
sm_perox_GSE71723v14
perox_vstc_GSE71723v14 <- vstc_dds_GSE71723v14[perox_GSE71723v14_de,]
dim(perox_vstc_GSE71723v14)
DEperox_zscore_GSE71723v14 <- (perox_vstc_GSE71723v14-rowMeans(perox_vstc_GSE71723v14))/(rowSds(as.matrix(perox_vstc_GSE71723v14)))
head(DEperox_zscore_GSE71723v14)
GSE71723v14_coln <- data.frame(sample = c("control","control","control","control","drought","drought","drought","drought"))
GSE71723v14_coln
row.names(GSE71723v14_coln) <- colnames(vstc_dds_GSE71723v14)
pheatmap(DEperox_zscore_GSE71723v14,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE71723v14_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE71723v14,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)

