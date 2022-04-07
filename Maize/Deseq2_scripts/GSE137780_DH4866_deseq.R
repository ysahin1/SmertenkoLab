library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
GSE137780_DH4866_raw_counts <- read.table("../GSE137780_DH4866_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE137780_DH4866_raw_counts)
row.names(GSE137780_DH4866_raw_counts) <- GSE137780_DH4866_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE137780_DH4866 <- GSE137780_DH4866_raw_counts$Length
names(gene_length_GSE137780_DH4866) <- GSE137780_DH4866_raw_counts$Geneid
##################################################
GSE137780_DH4866_raw_counts <- GSE137780_DH4866_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(GSE137780_DH4866_raw_counts) <- c("SRR10153104",
                        "SRR10153105",
                        "SRR10153116",
                        "SRR10153117")
GSE137780_DH4866_grow_CG <- c("control","control","drought","drought")
#GSE137780_DH4866_grow_CG <- factor(GSE137780_DH4866_grow_CG, levels = c("0","1","2"))
#GSE137780_DH4866_grow_CG <- relevel(GSE137780_DH4866_grow_CG, ref = "0")
####################################
### generate times series deseq object
####################################
GSE137780_DH4866_sample_info <- data.frame(condition=GSE137780_DH4866_grow_CG,
                          row.names=names(GSE137780_DH4866_raw_counts))

dds_GSE137780_DH4866 <- DESeqDataSetFromMatrix(countData = GSE137780_DH4866_raw_counts,
                                   colData = GSE137780_DH4866_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
#assay(DESeq.ds) %>% head

##################################
### filter out count datas and normalization
################################## 
dds_GSE137780_DH4866 <- dds_GSE137780_DH4866[rowSums(counts(dds_GSE137780_DH4866))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE137780_DH4866)) # number of genes counted
head(assay(dds_GSE137780_DH4866))

######normalization
dds_GSE137780_DH4866 <- estimateSizeFactors(dds_GSE137780_DH4866)
sizeFactors(dds_GSE137780_DH4866)
sfn_dds_GSE137780_DH4866 <- counts(dds_GSE137780_DH4866, normalized = TRUE)

########transformation
ln_dds_GSE137780_DH4866 <- log2(sfn_dds_GSE137780_DH4866+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE137780_DH4866, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE137780_DH4866 , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE137780_DH4866 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE137780_DH4866 <- vst(dds_GSE137780_DH4866, blind = TRUE)
vstc_dds_GSE137780_DH4866 <- assay(vst_dds_GSE137780_DH4866)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE137780_DH4866,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE137780_DH4866,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE137780_DH4866 <- as.dist(1- cor(vstc_dds_GSE137780_DH4866, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE137780_DH4866),
      labels = colnames(dis_vstc_dds_GSE137780_DH4866),
      xlab = "",
      main = "Hierarchical Clustering of GSE137780_DH4866")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE137780_DH4866))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE137780_DH4866)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE137780_DH4866)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE137780_DH4866")
print (P)

############## DEG analysis ############################
dds_GSE137780_DH4866 <- DESeq(dds_GSE137780_DH4866)
dds_GSE137780_DH4866 <- estimateSizeFactors(dds_GSE137780_DH4866)
dds_GSE137780_DH4866 <- estimateDispersions(dds_GSE137780_DH4866)
dds_GSE137780_DH4866 <- nbinomWaldTest(dds_GSE137780_DH4866)

dds_GSE137780_DH4866_results <- results(dds_GSE137780_DH4866, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE137780_DH4866_results)

dn_GSE137780_DH4866 <- rownames(subset(dds_GSE137780_DH4866_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE137780_DH4866)
dn_GSE137780_DH4866[c("Zm00001d017072", "Zm00001d003300")]
