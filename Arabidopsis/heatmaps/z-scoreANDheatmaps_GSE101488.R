library(NMF)
library(pheatmap)
peroxisome_gene_table <- read.delim("~/Desktop/MEGA/peroxisome-C3vsC4/DEG/GSE134945/peroxisome_geneID_sembol", header=FALSE)
peroxisome_gene_table <- as.data.frame(peroxisome_gene_table)
names(peroxisome_gene_table) <- c("geneSym", "gene_id")

perox_GSE101488_de <- dn_GSE101488[which(dn_GSE101488 %in% peroxisome_gene_table$gene_id)]
length(perox_GSE101488_de)


sm_perox_GSE101488 <- as.vector(peroxisome_gene_table[which(as.vector(peroxisome_gene_table$gene_id)%in% perox_GSE101488_de),1])
sm_perox_GSE101488
perox_vstc_GSE101488 <- vstc_dds_GSE101488[perox_GSE101488_de,]
DEperox_zscore_GSE101488 <- (perox_vstc_GSE101488-rowMeans(perox_vstc_GSE101488))/(rowSds(as.matrix(perox_vstc_GSE101488)))
head(DEperox_zscore_GSE101488)
GSE101488_coln <- data.frame(sample = c("control","control","control","ABA","ABA","ABA","control","control","control","ABA","ABA","ABA"))
GSE101488_coln
row.names(GSE101488_coln) <- colnames(vstc_dds_GSE101488)
pheatmap(DEperox_zscore_GSE101488,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE101488_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE101488,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)

