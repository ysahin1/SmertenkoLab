library(NMF)
library(pheatmap)
os_homo <- read.csv("../rice_peroxisome_homologs.csv", header = T, sep = "\t")
perox_GSE57950p_de <- dn_GSE57950p[which(dn_GSE57950p %in% os_homo$Os.locus)]
nperDEG <- length(perox_GSE57950p_de)
nperDEG
sym_idx <- match(perox_GSE57950p_de,os_homo$Os.locus)
sm_perox_GSE57950p <- as.vector(os_homo$Gene.name[sym_idx])
perox_vstc_GSE57950p <- vstc_dds_GSE57950p[perox_GSE57950p_de,]
dim(perox_vstc_GSE57950p)
DEperox_zscore_GSE57950p <- (perox_vstc_GSE57950p-rowMeans(perox_vstc_GSE57950p))/(rowSds(as.matrix(perox_vstc_GSE57950p)))
head(DEperox_zscore_GSE57950p)
GSE57950p_coln <- data.frame(sample = c("control", "control","drought","drought"))
GSE57950p_coln
row.names(GSE57950p_coln) <- colnames(vstc_dds_GSE57950p)
tiff("GSE57950hz_per.tiff", units="in", width=10, height=20, res=600)
pheatmap(DEperox_zscore_GSE57950p,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE57950p_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE57950p,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F, fontsize = 14)
dev.off()
