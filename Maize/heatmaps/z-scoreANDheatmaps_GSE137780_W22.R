library(pheatmap)
library(DESeq2)
zm_per <- read.delim("../os_zm_sym_per_homologs.txt")
perox_GSE137780_W22_de <- dn_GSE137780_W22[which(dn_GSE137780_W22 %in% zm_per$zm_locus)]
nperDEG <- length(perox_GSE137780_W22_de)
nperDEG
sym_idx <- match(perox_GSE137780_W22_de,zm_per$zm_locus)
sm_perox_GSE137780_W22 <- as.vector(zm_per$symbols[sym_idx])
sm_perox_GSE137780_W22
perox_vstc_GSE137780_W22 <- vstc_dds_GSE137780_W22[perox_GSE137780_W22_de,]
dim(perox_vstc_GSE137780_W22)
DEperox_zscore_GSE137780_W22 <- (perox_vstc_GSE137780_W22-rowMeans(perox_vstc_GSE137780_W22))/(rowSds(as.matrix(perox_vstc_GSE137780_W22)))
head(DEperox_zscore_GSE137780_W22)
GSE137780_W22_coln <- data.frame(sample = c("control","control","drought","drought"))
GSE137780_W22_coln
row.names(GSE137780_W22_coln) <- colnames(vstc_dds_GSE137780_W22)
tiff(filename = "GSE137780_W22.tiff", units = "in", width = 8, height = 16, res = 600)
pheatmap(DEperox_zscore_GSE137780_W22,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE137780_W22_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE137780_W22,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)
dev.off()
