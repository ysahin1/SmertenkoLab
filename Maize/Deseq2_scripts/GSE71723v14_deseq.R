library(magrittr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(DESeq2)


###################################################
## load file
###################################################
GSE71723v14_raw_counts <- read.table("./GSE71723v14_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE71723v14_raw_counts)
row.names(GSE71723v14_raw_counts) <- GSE71723v14_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE71723v14 <- GSE71723v14_raw_counts$Length
names(gene_length_GSE71723v14) <- GSE71723v14_raw_counts$Geneid
##################################################
GSE71723v14_raw_counts <- GSE71723v14_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(GSE71723v14_raw_counts) <- c("GSM1843748",
                        "GSM1843749",
                        "GSM1843750",
                        "GSM1843751",
                        "GSM1843780",
                        "GSM1843781",
                        "GSM1843782",
                        "GSM1843783")
GSE71723v14_grow_CG <- c("control","control","control","control","drought","drought","drought","drought")
#GSE71723v14_grow_CG <- factor(GSE71723v14_grow_CG, levels = c("control","drought"))
#GSE71723v14_grow_CG <- relevel(GSE71723v14_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
GSE71723v14_sample_info <- data.frame(condition=GSE71723v14_grow_CG,
                          row.names=names(GSE71723v14_raw_counts))

dds_GSE71723v14 <- DESeqDataSetFromMatrix(countData = GSE71723v14_raw_counts,
                                   colData = GSE71723v14_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
#assay(DESeq.ds) %>% head

##################################
### filter out count datas and normalization
################################## 
dds_GSE71723v14 <- dds_GSE71723v14[rowSums(counts(dds_GSE71723v14))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE71723v14)) # number of genes counted
head(assay(dds_GSE71723v14))

######normalization
dds_GSE71723v14 <- estimateSizeFactors(dds_GSE71723v14)
sizeFactors(dds_GSE71723v14)
sfn_dds_GSE71723v14 <- counts(dds_GSE71723v14, normalized = TRUE)

########transformation
ln_dds_GSE71723v14 <- log2(sfn_dds_GSE71723v14+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE71723v14, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE71723v14 , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE71723v14 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE71723v14 <- vst(dds_GSE71723v14, blind = TRUE)
vstc_dds_GSE71723v14 <- assay(vst_dds_GSE71723v14)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE71723v14,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE71723v14,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE71723v14 <- as.dist(1- cor(vstc_dds_GSE71723v14, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE71723v14),
      labels = colnames(dis_vstc_dds_GSE71723v14),
      xlab = "",
      main = "Hierarchical Clustering of GSE71723v14")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE71723v14))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE71723v14)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE71723v14)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE71723v14")
print (P)

############## DEG analysis ############################
dds_GSE71723v14 <- DESeq(dds_GSE71723v14)
dds_GSE71723v14 <- estimateSizeFactors(dds_GSE71723v14)
dds_GSE71723v14 <- estimateDispersions(dds_GSE71723v14)
dds_GSE71723v14 <- nbinomWaldTest(dds_GSE71723v14)

dds_GSE71723v14_results <- results(dds_GSE71723v14, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE71723v14_results)

dn_GSE71723v14 <- rownames(subset(dds_GSE71723v14_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE71723v14)
