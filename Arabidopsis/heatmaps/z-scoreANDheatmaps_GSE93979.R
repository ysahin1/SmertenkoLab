library(NMF)
library(pheatmap)
peroxisome_gene_table <- read.delim("~/Desktop/MEGA/peroxisome-C3vsC4/DEG/GSE134945/peroxisome_geneID_sembol", header=FALSE)
peroxisome_gene_table <- as.data.frame(peroxisome_gene_table)
names(peroxisome_gene_table) <- c("geneSym", "gene_id")

perox_GSE93979_de <- dn_GSE93979[which(dn_GSE93979 %in% peroxisome_gene_table$gene_id)]
length(perox_GSE93979_de)


sm_perox_GSE93979 <- as.vector(peroxisome_gene_table[which(as.vector(peroxisome_gene_table$gene_id)%in% perox_GSE93979_de),1])
sm_perox_GSE93979
perox_vstc_GSE93979 <- vstc_dds_GSE93979[perox_GSE93979_de,]
DEperox_zscore_GSE93979 <- (perox_vstc_GSE93979-rowMeans(perox_vstc_GSE93979))/(rowSds(as.matrix(perox_vstc_GSE93979)))
head(DEperox_zscore_GSE93979)
GSE93979_coln <- data.frame(sample = rep(c("control", "drought"), c(2,2)))
GSE93979_coln
row.names(GSE93979_coln) <- colnames(vstc_dds_GSE93979)
pheatmap(DEperox_zscore_GSE93979,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE93979_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE93979,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)

