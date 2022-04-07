library(NMF)
library(pheatmap)
zm_per <- read.delim("../os_zm_sym_per_homologs.txt")
perox_SRP132035_de <- dn_SRP132035[which(dn_SRP132035 %in% zm_per$zm_locus)]
nperDEG <- length(perox_SRP132035_de)
nperDEG
sym_idx <- match(perox_SRP132035_de,zm_per$zm_locus)
sm_perox_SRP132035 <- as.vector(zm_per$symbols[sym_idx])
sm_perox_SRP132035
perox_vstc_SRP132035 <- vstc_dds_SRP132035[perox_SRP132035_de,]
dim(perox_vstc_SRP132035)
DEperox_zscore_SRP132035 <- (perox_vstc_SRP132035-rowMeans(perox_vstc_SRP132035))/(rowSds(as.matrix(perox_vstc_SRP132035)))
head(DEperox_zscore_SRP132035)
SRP132035_coln <- data.frame(sample = c("control","control",
                                        "drought","drought","drought","drought","drought","drought",
                                        "control","control","control","control"))
SRP132035_coln
row.names(SRP132035_coln) <- colnames(vstc_dds_SRP132035)
tiff(filename = "SRP132035.tiff", units = "in", width = 8, height = 16, res = 600)
pheatmap(DEperox_zscore_SRP132035,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = SRP132035_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_SRP132035,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)
dev.off()
