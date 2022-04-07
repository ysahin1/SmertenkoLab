library(NMF)
library(pheatmap)
library(dplyr)
peroxisome_gene_table <- read.delim("~/Desktop/MEGA/peroxisome-C3vsC4/DEG/GSE134945/peroxisome_geneID_sembol", header=FALSE)
peroxisome_gene_table <- as.data.frame(peroxisome_gene_table)
names(peroxisome_gene_table) <- c("geneSym", "gene_id")
row.names(peroxisome_gene_table) <- peroxisome_gene_table$gene_id

perox_GSE134945_de <- dn_GSE134945[which(dn_GSE134945 %in% peroxisome_gene_table$gene_id)]
perox_GSE134945_de
as.vector(peroxisome_gene_table$gene_id)


sm_perox_GSE134945 <- as.vector(peroxisome_gene_table[which(as.vector(peroxisome_gene_table$gene_id)%in% perox_GSE134945_de),1])
sm_perox_GSE134945
head(vstc_dds_GSE134945)
perox_vstc_GSE134945 <- vstc_dds_GSE134945[perox_GSE134945_de,]
DEperox_zscore_GSE134945 <- (perox_vstc_GSE134945-rowMeans(perox_vstc_GSE134945))/(rowSds(as.matrix(perox_vstc_GSE134945)))
head(DEperox_zscore_GSE134945)
png("S1.png", width=800, height=600)
pheatmap(DEperox_zscore_GSE134945,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = my_sample_col,
         cex=1, clustering_distance_rows="euclidean", labels_row =  sm_perox_GSE134945,
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)
dev.off()
#### bütün peroxizome genleri için heatmap DEGler dahil
dim(vstc_dds_GSE134945)
all_perox_vstc_GSE134945 <- vstc_dds_GSE134945[which(rownames(vstc_dds_GSE134945) %in% perox_gene_id),]
all_perox_zscore_GSE134945 <- (all_perox_vstc_GSE134945-rowMeans(all_perox_vstc_GSE134945))/(rowSds(as.matrix(all_perox_vstc_GSE134945)))
png("S1.png", width=800, height=1200)
#par(mar=c(5,5,9,24))
par(cex=1, mai=c(0.5,0.5,0.5,0.5))
pheatmap(all_perox_vstc_GSE134945,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",
         clustering_distance_rows="euclidean", 
         fontsize = 18, 
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F, 
         labels_row = as.character(rownames(all_perox_vstc_GSE134945)),
         annotation_col = my_sample_col)
dev.off()
my_sample_col <- data.frame(sample = rep(c("control", "drought"), c(3,3)))
my_sample_col
row.names(my_sample_col) <- colnames(all_perox_vstc_GSE134945)

