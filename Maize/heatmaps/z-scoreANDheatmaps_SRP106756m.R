library(pheatmap)
zm_per <- read.delim("../os_zm_sym_per_homologs.txt")
perox_SRP106756_de <- dn_SRP106756[which(dn_SRP106756 %in% zm_per$zm_locus)]
nperDEG <- length(perox_SRP106756_de)
nperDEG
sym_idx <- match(perox_SRP106756_de,zm_per$zm_locus)
sm_perox_SRP106756 <- as.vector(zm_per$symbols[sym_idx])
sm_perox_SRP106756
perox_vstc_SRP106756 <- vstc_dds_SRP106756[perox_SRP106756_de,]
dim(perox_vstc_SRP106756)
DEperox_zscore_SRP106756 <- (perox_vstc_SRP106756-rowMeans(perox_vstc_SRP106756))/(rowSds(as.matrix(perox_vstc_SRP106756)))
head(DEperox_zscore_SRP106756)
SRP106756_coln <- data.frame(sample = c("0","1","2","3","0","1","2","3"))
SRP106756_coln
row.names(SRP106756_coln) <- colnames(vstc_dds_SRP106756)
tiff(filename = "SRP106756.tiff", units = "in", width = 8, height = 16, res=600)
pheatmap(DEperox_zscore_SRP106756,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = SRP106756_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_SRP106756,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)
dev.off()
