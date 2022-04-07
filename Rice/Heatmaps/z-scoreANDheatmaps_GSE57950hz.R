library(NMF)
library(pheatmap)
os_homo <- read.csv("../rice_peroxisome_homologs.csv", header = T, sep = "\t")
perox_GSE57950hz_de <- dn_GSE57950hz[which(dn_GSE57950hz %in% os_homo$Os.locus)]
nperDEG <- length(perox_GSE57950hz_de)
nperDEG
sym_idx <- match(perox_GSE57950hz_de,os_homo$Os.locus)
sm_perox_GSE57950hz <- as.vector(os_homo$Gene.name[sym_idx])
perox_vstc_GSE57950hz <- vstc_dds_GSE57950hz[perox_GSE57950hz_de,]
dim(perox_vstc_GSE57950hz)
DEperox_zscore_GSE57950hz <- (perox_vstc_GSE57950hz-rowMeans(perox_vstc_GSE57950hz))/(rowSds(as.matrix(perox_vstc_GSE57950hz)))
head(DEperox_zscore_GSE57950hz)
GSE57950hz_coln <- data.frame(sample = c("control", "control","drought","drought"))
GSE57950hz_coln
row.names(GSE57950hz_coln) <- colnames(vstc_dds_GSE57950hz)
tiff("GSE57950hz_per.tiff", units="in", width=10, height=20, res=600)
pheatmap(DEperox_zscore_GSE57950hz,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE57950hz_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE57950hz,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F, fontsize = 14)
dev.off()

