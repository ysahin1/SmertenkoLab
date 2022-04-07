library(genefilter)
library(pheatmap)
autophagy_gene_table <- read.delim("/media/yunus/TOSHIBA1TB/c3-c4/peroxisome-C3vsC4/autophagy/autophagy_homologs.csv", header=T, stringsAsFactors = F)
dim(autophagy_gene_table)
locus_IDs <- strsplit(as.vector(autophagy_gene_table$solanum_tuberosum), "," )
nIDs <- sapply(1:length(locus_IDs), function(x) length(locus_IDs[[x]]))
locus_IDs <- gsub(" ", "", unlist(locus_IDs), fixed = TRUE)
os_autophagy_genes <- data.frame(gene_symbol = rep(autophagy_gene_table$gene_symbol, nIDs), 
                                 locus_IDs = unlist(locus_IDs), 
                                 stringsAsFactors=FALSE )
os_autophagy_genes$locus_IDs <- gsub("DMP", "DMG", os_autophagy_genes$locus_IDs)
os_autophagy_genes <- os_autophagy_genes[!duplicated(os_autophagy_genes$locus_IDs), ]
#####list comes from pythozome, no need to filter duplicates.
#autophagy_gene_table <- autophagy_gene_table[!duplicated(autophagy_gene_table$at_locus_id),]
#autophagy_gene_table$at_locus_id <- toupper(autophagy_gene_table$at_locus_id)

autophagy_GSE97776_tajfun_de <- dn_GSE97776_tajfun[which(dn_GSE97776_tajfun %in% os_autophagy_genes$locus_IDs)]
length(autophagy_GSE97776_tajfun_de)

sm_atg_GSE97776_tajfun <- as.vector(os_autophagy_genes[which(as.vector(os_autophagy_genes$locus_IDs)%in% autophagy_GSE97776_tajfun_de),])
dim(sm_atg_GSE97776_tajfun)

autophagy_vstc_GSE97776_tajfun <- vstc_dds_GSE97776_tajfun[autophagy_GSE97776_tajfun_de,]
DEauto_zscore_GSE97776_tajfun <- (autophagy_vstc_GSE97776_tajfun-rowMeans(autophagy_vstc_GSE97776_tajfun))/(rowSds(as.matrix(autophagy_vstc_GSE97776_tajfun)))
head(DEauto_zscore_GSE97776_tajfun)

GSE97776_tajfun_coln <- data.frame(sample = GSE97776_tajfun_grow_CG)
GSE97776_tajfun_coln
row.names(GSE97776_tajfun_coln) <- colnames(vstc_dds_GSE97776_tajfun)

ind_idx <- match(rownames(DEauto_zscore_GSE97776_tajfun),sm_atg_GSE97776_tajfun$locus_IDs)
tiff("GSE97776_tajfun_auto.tiff", units="in", width=8, height=10, res=600)
pheatmap(DEauto_zscore_GSE97776_tajfun,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE97776_tajfun_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  paste(sm_atg_GSE97776_tajfun$gene_symbol[ind_idx],sm_atg_GSE97776_tajfun$locus_IDs[ind_idx] ,sep = "-"),
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)
dev.off()

