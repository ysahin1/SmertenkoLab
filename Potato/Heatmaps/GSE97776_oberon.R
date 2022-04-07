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

autophagy_GSE97776_oberon_de <- dn_GSE97776_oberon[which(dn_GSE97776_oberon %in% os_autophagy_genes$locus_IDs)]
length(autophagy_GSE97776_oberon_de)

sm_atg_GSE97776_oberon <- as.vector(os_autophagy_genes[which(as.vector(os_autophagy_genes$locus_IDs)%in% autophagy_GSE97776_oberon_de),])
dim(sm_atg_GSE97776_oberon)

autophagy_vstc_GSE97776_oberon <- vstc_dds_GSE97776_oberon[autophagy_GSE97776_oberon_de,]
DEauto_zscore_GSE97776_oberon <- (autophagy_vstc_GSE97776_oberon-rowMeans(autophagy_vstc_GSE97776_oberon))/(rowSds(as.matrix(autophagy_vstc_GSE97776_oberon)))
head(DEauto_zscore_GSE97776_oberon)

GSE97776_oberon_coln <- data.frame(sample = GSE97776_oberon_grow_CG)
GSE97776_oberon_coln
row.names(GSE97776_oberon_coln) <- colnames(vstc_dds_GSE97776_oberon)

ind_idx <- match(rownames(DEauto_zscore_GSE97776_oberon),sm_atg_GSE97776_oberon$locus_IDs)
tiff("GSE97776_oberon_auto.tiff", units="in", width=8, height=10, res=600)
pheatmap(DEauto_zscore_GSE97776_oberon,
         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",annotation_col = GSE97776_oberon_coln,
         cex=1, clustering_distance_rows="euclidean", labels_row =  paste(sm_atg_GSE97776_oberon$gene_symbol[ind_idx],sm_atg_GSE97776_oberon$locus_IDs[ind_idx] ,sep = "-"),
         clustering_distance_cols="euclidean", clustering_method="complete", border_color=F)
dev.off()

