library(pheatmap)
library(DESeq2)
zm_per <- read.delim("../os_zm_sym_per_homologs.txt")
perox_GSE48507_de <- dn_GSE48507[which(dn_GSE48507 %in% zm_per$zm_locus)]
nperDEG <- length(perox_GSE48507_de)
nperDEG
sym_idx <- match(perox_GSE48507_de,zm_per$zm_locus)
sm_perox_GSE48507 <- as.vector(zm_per$symbols[sym_idx])
sm_perox_GSE48507
perox_vstc_GSE48507 <- vstc_dds_GSE48507[perox_GSE48507_de,]
dim(perox_vstc_GSE48507)
DEperox_zscore_GSE48507 <- (perox_vstc_GSE48507-rowMeans(perox_vstc_GSE48507))/(rowSds(as.matrix(perox_vstc_GSE48507)))
head(DEperox_zscore_GSE48507)
GSE48507_coln <- data.frame(sample = c("control","control","drought_s1","drought_s1","drought_s2","drought_s2"))
GSE48507_coln
row.names(GSE48507_coln) <- colnames(vstc_dds_GSE48507)
tiff(filename = "GSE48507.tiff", units = "in", width = 10, height = 16, res=600)
pheatmap(DEperox_zscore_GSE48507,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE48507_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE48507,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F, fontsize = 14)
dev.off()
