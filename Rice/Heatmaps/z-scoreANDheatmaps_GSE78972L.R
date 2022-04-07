library(NMF)
library(pheatmap)
os_homo <- read.csv("../rice_peroxisome_homologs.csv", header = T, sep = "\t")
perox_GSE78972l_de <- dn_GSE78972l[which(dn_GSE78972l %in% os_homo$Os.locus)]
nperDEG <- length(perox_GSE78972l_de)
nperDEG
sym_idx <- match(perox_GSE78972l_de,os_homo$Os.locus)
sm_perox_GSE78972l <- as.vector(os_homo$Gene.name[sym_idx])
perox_vstc_GSE78972l <- vstc_dds_GSE78972l[perox_GSE78972l_de,]
dim(perox_vstc_GSE78972l)
DEperox_zscore_GSE78972l <- (perox_vstc_GSE78972l-rowMeans(perox_vstc_GSE78972l))/(rowSds(as.matrix(perox_vstc_GSE78972l)))
head(DEperox_zscore_GSE78972l)
GSE78972l_coln <- data.frame(sample = c("control", "control","drought","drought","control","drought"))
GSE78972l_coln
row.names(GSE78972l_coln) <- colnames(vstc_dds_GSE78972l)
tiff("GSE78972l_per.tiff", units="in", width=10, height=25, res=600)
pheatmap(DEperox_zscore_GSE78972l,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE78972l_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE78972l,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F, fontsize = 14)

dev.off()