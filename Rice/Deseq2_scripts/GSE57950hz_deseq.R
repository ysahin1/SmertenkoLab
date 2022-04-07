library(ggplot2)
library(DESeq2)

###################################################
## load file
###################################################
GSE57950hz_raw_counts <- read.table("raw_counts/HHZ_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE57950hz_raw_counts)
row.names(GSE57950hz_raw_counts) <- GSE57950hz_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE57950hz <- GSE57950hz_raw_counts$Length
names(gene_length_GSE57950hz) <- GSE57950hz_raw_counts$Geneid
##################################################
GSE57950hz_raw_counts <- GSE57950hz_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(GSE57950hz_raw_counts) <- c("GSM1398429",
                        "GSM1398430",
                        "GSM1398441",
                        "GSM1398442")
GSE57950hz_grow_CG <- c("control","control","drought","drought")
### generate times series deseq object
####################################
GSE57950hz_sample_info <- data.frame(condition=GSE57950hz_grow_CG,
                          row.names=names(GSE57950hz_raw_counts))

dds_GSE57950hz <- DESeqDataSetFromMatrix(countData = GSE57950hz_raw_counts,
                                   colData = GSE57950hz_sample_info,
                                   design = ~ condition)
##################################
### filter out count datas and normalization
################################## 
dds_GSE57950hz <- dds_GSE57950hz[rowSums(counts(dds_GSE57950hz))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE57950hz)) # number of genes counted
head(assay(dds_GSE57950hz))

######normalization
dds_GSE57950hz <- estimateSizeFactors(dds_GSE57950hz)
sizeFactors(dds_GSE57950hz)
sfn_dds_GSE57950hz <- counts(dds_GSE57950hz, normalized = TRUE)

########transformation
ln_dds_GSE57950hz <- log2(sfn_dds_GSE57950hz+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE57950hz, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE57950hz , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE57950hz [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE57950hz <- vst(dds_GSE57950hz, blind = TRUE)
vstc_dds_GSE57950hz <- assay(vst_dds_GSE57950hz)

######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE57950hz,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE57950hz,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE57950hz <- as.dist(1- cor(vstc_dds_GSE57950hz, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE57950hz),
      labels = colnames(dis_vstc_dds_GSE57950hz),
      xlab = "",
      main = "Hierarchical Clustering of GSE57950hz")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE57950hz))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE57950hz)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE57950hz)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE57950hz")
print (P)

############## DEG analysis ############################
dds_GSE57950hz <- DESeq(dds_GSE57950hz)

dds_GSE57950hz_results <- results(dds_GSE57950hz, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE57950hz_results)

dn_GSE57950hz <- rownames(subset(dds_GSE57950hz_results, 
                             padj < 0.05 & abs(log2FoldChange) > 2))
length(dn_GSE57950hz)
