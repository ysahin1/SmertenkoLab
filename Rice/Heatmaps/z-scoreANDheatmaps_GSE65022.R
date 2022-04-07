library(NMF)
library(pheatmap)
os_homo <- read.csv("../rice_peroxisome_homologs.csv", header = T, sep = "\t")
perox_GSE65022_de <- dn_GSE65022[which(dn_GSE65022 %in% os_homo$Os.locus)]
nperDEG <- length(perox_GSE65022_de)
nperDEG
sym_idx <- match(perox_GSE65022_de,os_homo$Os.locus)
sm_perox_GSE65022 <- as.vector(os_homo$Gene.name[sym_idx])
perox_vstc_GSE65022 <- vstc_dds_GSE65022[perox_GSE65022_de,]
dim(perox_vstc_GSE65022)
DEperox_zscore_GSE65022 <- (perox_vstc_GSE65022-rowMeans(perox_vstc_GSE65022))/(rowSds(as.matrix(perox_vstc_GSE65022)))
head(DEperox_zscore_GSE65022)
GSE65022_coln <- data.frame(sample = c("drought", "drought","control","control"))
GSE65022_coln
row.names(GSE65022_coln) <- colnames(vstc_dds_GSE65022)
tiff("GSE65022_per.tiff", units="in", width=10, height=20, res=600)
pheatmap(DEperox_zscore_GSE65022,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE65022_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE65022,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F, fontsize = 14)
dev.off()
