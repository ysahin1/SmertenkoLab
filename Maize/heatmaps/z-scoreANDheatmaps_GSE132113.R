library(NMF)
library(pheatmap)
zm_per <- read.delim("os_zm_sym_per_homologs.txt")
perox_GSE132113_de <- dn_GSE132113[which(dn_GSE132113 %in% zm_per$zm_locus)]
nperDEG <- length(perox_GSE132113_de)
nperDEG
sym_idx <- match(perox_GSE132113_de,zm_per$zm_locus)
sm_perox_GSE132113 <- as.vector(zm_per$symbols[sym_idx])
sm_perox_GSE132113
perox_vstc_GSE132113 <- vstc_dds_GSE132113[perox_GSE132113_de,]
dim(perox_vstc_GSE132113)
DEperox_zscore_GSE132113 <- (perox_vstc_GSE132113-rowMeans(perox_vstc_GSE132113))/(rowSds(as.matrix(perox_vstc_GSE132113)))
head(DEperox_zscore_GSE132113)
GSE132113_coln <- data.frame(sample = c("control", "control","drought","drought"))
GSE132113_coln
row.names(GSE132113_coln) <- colnames(vstc_dds_GSE132113)
pheatmap(DEperox_zscore_GSE132113,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE132113_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE132113,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)

