library(NMF)
library(pheatmap)
zm_per <- read.delim("os_zm_sym_per_homologs.txt")
perox_GSE40070_de <- dn_GSE40070[which(dn_GSE40070 %in% zm_per$zm_locus)]
nperDEG <- length(perox_GSE40070_de)
nperDEG
sym_idx <- match(perox_GSE40070_de,zm_per$zm_locus)
sm_perox_GSE40070 <- as.vector(zm_per$symbols[sym_idx])
sm_perox_GSE40070
perox_vstc_GSE40070 <- vstc_dds_GSE40070[perox_GSE40070_de,]
dim(perox_vstc_GSE40070)
DEperox_zscore_GSE40070 <- (perox_vstc_GSE40070-rowMeans(perox_vstc_GSE40070))/(rowSds(as.matrix(perox_vstc_GSE40070)))
head(DEperox_zscore_GSE40070)
GSE40070_coln <- data.frame(sample = c("control", "control","drought","drought"))
GSE40070_coln
row.names(GSE40070_coln) <- colnames(vstc_dds_GSE40070)
pheatmap(DEperox_zscore_GSE40070,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE40070_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE40070,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)

