library(magrittr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(DESeq2)


###################################################
## load file
###################################################
GSE71723_raw_counts <- read.table("./GSE71723_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE71723_raw_counts)
row.names(GSE71723_raw_counts) <- GSE71723_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE71723 <- GSE71723_raw_counts$Length
names(gene_length_GSE71723) <- GSE71723_raw_counts$Geneid
##################################################
GSE71723_raw_counts <- GSE71723_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(GSE71723_raw_counts) <- c("GSM1843740",
                        "GSM1843741",
                        "GSM1843742",
                        "GSM1843743",
                        "GSM1843772",
                        "GSM1843773",
                        "GSM1843774",
                        "GSM1843775")
GSE71723_grow_CG <- c("control","control","control","control","drought","drought","drought","drought")
#GSE71723_grow_CG <- factor(GSE71723_grow_CG, levels = c("control","drought"))
#GSE71723_grow_CG <- relevel(GSE71723_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
GSE71723_sample_info <- data.frame(condition=GSE71723_grow_CG,
                          row.names=names(GSE71723_raw_counts))

dds_GSE71723 <- DESeqDataSetFromMatrix(countData = GSE71723_raw_counts,
                                   colData = GSE71723_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
#assay(DESeq.ds) %>% head

##################################
### filter out count datas and normalization
################################## 
dds_GSE71723 <- dds_GSE71723[rowSums(counts(dds_GSE71723))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE71723)) # number of genes counted
head(assay(dds_GSE71723))

######normalization
dds_GSE71723 <- estimateSizeFactors(dds_GSE71723)
sizeFactors(dds_GSE71723)
sfn_dds_GSE71723 <- counts(dds_GSE71723, normalized = TRUE)

########transformation
ln_dds_GSE71723 <- log2(sfn_dds_GSE71723+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE71723, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE71723 , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE71723 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE71723 <- vst(dds_GSE71723, blind = TRUE)
vstc_dds_GSE71723 <- assay(vst_dds_GSE71723)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE71723,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE71723,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE71723 <- as.dist(1- cor(vstc_dds_GSE71723, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE71723),
      labels = colnames(dis_vstc_dds_GSE71723),
      xlab = "",
      main = "Hierarchical Clustering of GSE71723")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE71723))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE71723)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE71723)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE71723")
print (P)

############## DEG analysis ############################
dds_GSE71723 <- DESeq(dds_GSE71723)
dds_GSE71723 <- estimateSizeFactors(dds_GSE71723)
dds_GSE71723 <- estimateDispersions(dds_GSE71723)
dds_GSE71723 <- nbinomWaldTest(dds_GSE71723)

dds_GSE71723_results <- results(dds_GSE71723, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE71723_results)

dn_GSE71723 <- rownames(subset(dds_GSE71723_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE71723)
