library(NMF)
library(pheatmap)
library(DESeq2)
os_homo <- read.csv("../rice_peroxisome_homologs.csv", header = T, sep = "\t")
perox_GSE57950h_de <- dn_GSE57950h[which(dn_GSE57950h %in% os_homo$Os.locus)]
nperDEG <- length(perox_GSE57950h_de)
nperDEG
sym_idx <- match(perox_GSE57950h_de,os_homo$Os.locus)
sm_perox_GSE57950h <- as.vector(os_homo$Gene.name[sym_idx])
perox_vstc_GSE57950h <- vstc_dds_GSE57950h[perox_GSE57950h_de,]
dim(perox_vstc_GSE57950h)
DEperox_zscore_GSE57950h <- (perox_vstc_GSE57950h-rowMeans(perox_vstc_GSE57950h))/(rowSds(as.matrix(perox_vstc_GSE57950h)))
head(DEperox_zscore_GSE57950h)
GSE57950h_coln <- data.frame(sample = c("control", "control","drought","drought"))
GSE57950h_coln
row.names(GSE57950h_coln) <- colnames(vstc_dds_GSE57950h)
tiff("GSE57950h.tiff", units="in", width=10, height=20, res=600)
pheatmap(DEperox_zscore_GSE57950h,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE57950h_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE57950h,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F, fontsize = 14)
dev.off()
  