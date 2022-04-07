library(genefilter)
library(pheatmap)
autophagy_gene_table <- read.delim("/media/yunus/TOSHIBA1TB/c3-c4/peroxisome-C3vsC4/autophagy/autophagy_homologs.csv", header=T, stringsAsFactors = F)
dim(autophagy_gene_table)
locus_IDs <- strsplit(as.vector(autophagy_gene_table$solanum_lycopersicum), "," )
nIDs <- sapply(1:length(locus_IDs), function(x) length(locus_IDs[[x]]))
locus_IDs <- gsub(" ", "", unlist(locus_IDs), fixed = TRUE)
os_autophagy_genes <- data.frame(gene_symbol = rep(autophagy_gene_table$gene_symbol, nIDs), 
                                 locus_IDs = unlist(locus_IDs), 
                                 stringsAsFactors=FALSE )
#####list comes from pythozome, no need to filter duplicates.
#autophagy_gene_table <- autophagy_gene_table[!duplicated(autophagy_gene_table$at_locus_id),]
#autophagy_gene_table$at_locus_id <- toupper(autophagy_gene_table$at_locus_id)

autophagy_GSE151277_d3_de <- dn_GSE151277_d3[which(dn_GSE151277_d3 %in% os_autophagy_genes$locus_IDs)]
length(autophagy_GSE151277_d3_de)

sm_atg_GSE151277_d3 <- as.vector(os_autophagy_genes[which(as.vector(os_autophagy_genes$locus_IDs)%in% autophagy_GSE151277_d3_de),])
dim(sm_atg_GSE151277_d3)

autophagy_vstc_GSE151277_d3 <- vstc_dds_GSE151277_d3[autophagy_GSE151277_d3_de,]
DEauto_zscore_GSE151277_d3 <- (autophagy_vstc_GSE151277_d3-rowMeans(autophagy_vstc_GSE151277_d3))/(rowSds(as.matrix(autophagy_vstc_GSE151277_d3)))
head(DEauto_zscore_GSE151277_d3)

GSE151277_d3_coln <- data.frame(sample = GSE151277_d3_grow_CG)
GSE151277_d3_coln
row.names(GSE151277_d3_coln) <- colnames(vstc_dds_GSE151277_d3)

ind_idx <- match(rownames(DEauto_zscore_GSE151277_d3),sm_atg_GSE151277_d3$locus_IDs)
tiff("GSE151277_d3_auto.tiff", units="in", width=8, height=10, res=600)
pheatmap(DEauto_zscore_GSE151277_d3,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE151277_d3_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  paste(sm_atg_GSE151277_d3$gene_symbol[ind_idx],sm_atg_GSE151277_d3$locus_IDs[ind_idx] ,sep = "-"),
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)
dev.off()

