library(NMF)
library(pheatmap)
peroxisome_gene_table <- read.delim("~/Desktop/MEGA/peroxisome-C3vsC4/DEG/GSE134945/peroxisome_geneID_sembol", header=FALSE)
peroxisome_gene_table <- as.data.frame(peroxisome_gene_table)
names(peroxisome_gene_table) <- c("geneSym", "gene_id")

perox_GSE108610_de <- dn_GSE108610[which(dn_GSE108610 %in% peroxisome_gene_table$gene_id)]
length(perox_GSE108610_de)


sm_perox_GSE108610 <- as.vector(peroxisome_gene_table[which(as.vector(peroxisome_gene_table$gene_id)%in% perox_GSE108610_de),1])
sm_perox_GSE108610
perox_vstc_GSE108610 <- vstc_dds_GSE108610[perox_GSE108610_de,]
DEperox_zscore_GSE108610 <- (perox_vstc_GSE108610-rowMeans(perox_vstc_GSE108610))/(rowSds(as.matrix(perox_vstc_GSE108610)))
head(DEperox_zscore_GSE108610)
GSE108610_coln <- data.frame(sample = rep(c("control", "drought"), c(2,2)))

row.names(GSE108610_coln) <- colnames(vstc_dds_GSE108610)
pheatmap(DEperox_zscore_GSE108610,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE108610_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE108610,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)

