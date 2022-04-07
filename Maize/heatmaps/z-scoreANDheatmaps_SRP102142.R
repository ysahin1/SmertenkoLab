library(NMF)
library(pheatmap)
zm_per <- read.delim("os_zm_sym_per_homologs.txt")
perox_SRP102142_de <- dn_SRP102142[which(dn_SRP102142 %in% zm_per$zm_locus)]
nperDEG <- length(perox_SRP102142_de)
nperDEG
sym_idx <- match(perox_SRP102142_de,zm_per$zm_locus)
sm_perox_SRP102142 <- as.vector(zm_per$symbols[sym_idx])
sm_perox_SRP102142
perox_vstc_SRP102142 <- vstc_dds_SRP102142[perox_SRP102142_de,]
dim(perox_vstc_SRP102142)
DEperox_zscore_SRP102142 <- (perox_vstc_SRP102142-rowMeans(perox_vstc_SRP102142))/(rowSds(as.matrix(perox_vstc_SRP102142)))
head(DEperox_zscore_SRP102142)
SRP102142_coln <- data.frame(sample = c("control", "control","drought","drought"))
SRP102142_coln
row.names(SRP102142_coln) <- colnames(vstc_dds_SRP102142)
pheatmap(DEperox_zscore_SRP102142,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = SRP102142_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_SRP102142,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)

