library(pheatmap)
zm_per <- read.delim("../os_zm_sym_per_homologs.txt")
perox_GSE137780_b104_de <- dn_GSE137780_b104[which(dn_GSE137780_b104 %in% zm_per$zm_locus)]
nperDEG <- length(perox_GSE137780_b104_de)
nperDEG
sym_idx <- match(perox_GSE137780_b104_de,zm_per$zm_locus)
sm_perox_GSE137780_b104 <- as.vector(zm_per$symbols[sym_idx])
sm_perox_GSE137780_b104
perox_vstc_GSE137780_b104 <- vstc_dds_GSE137780_b104[perox_GSE137780_b104_de,]
dim(perox_vstc_GSE137780_b104)
DEperox_zscore_GSE137780_b104 <- (perox_vstc_GSE137780_b104-rowMeans(perox_vstc_GSE137780_b104))/(rowSds(as.matrix(perox_vstc_GSE137780_b104)))
head(DEperox_zscore_GSE137780_b104)
GSE137780_b104_coln <- data.frame(sample = c("control","control","drought","drought"))
GSE137780_b104_coln
row.names(GSE137780_b104_coln) <- colnames(vstc_dds_GSE137780_b104)
tiff(filename = "GSE137780_b104.tiff", units = "in", width = 8, height = 16, res = 600)
pheatmap(DEperox_zscore_GSE137780_b104,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE137780_b104_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE137780_b104,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)

dev.off()