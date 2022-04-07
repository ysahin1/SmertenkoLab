library(pheatmap)
sb_per <- read.delim("os_bsv3_peroxisome_homologs_edit.txt")

# HOW MANY PEROXİSOME GENES ARE DEGS
perox_GSE80699_de <- dn_GSE80699[which(dn_GSE80699 %in% as.vector(sb_per$bs_id))]
head(dn_GSE80699)
head(as.vector(sb_per$bs_id))
nperDEG <- length(perox_GSE80699_de)
nperDEG # IT GIVES NUMBER OF PEROXISOME DEGS
####### RETRIEVE SYMBOLS OF DEGS

sym_idx <- match(perox_GSE80699_de, sb_per$bs_id)
sm_perox_GSE80699 <- as.vector(sb_per$symbols[sym_idx])
sm_perox_GSE80699
perox_vstc_GSE80699 <- vstc_dds_GSE80699[perox_GSE80699_de,]
dim(perox_vstc_GSE80699)
DEperox_zscore_GSE80699 <- (perox_vstc_GSE80699-rowMeans(perox_vstc_GSE80699))/(rowSds(as.matrix(perox_vstc_GSE80699)))
head(DEperox_zscore_GSE80699)
GSE80699_coln <- data.frame(sample = c("control","control","control","drought","drought","drought"))
GSE80699_coln
row.names(GSE80699_coln) <- colnames(vstc_dds_GSE80699)
tiff(filename = "GSE80699_w3.tiff", units = "in", width = 8, height = 16, res=600)
pheatmap(DEperox_zscore_GSE80699,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE80699_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE80699,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F,fontsize = 14)
dev.off()