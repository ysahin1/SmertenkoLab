library(magrittr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(DESeq2)


###################################################
## load file
###################################################
GSE71046_raw_counts <- read.table("./GSE71046_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE71046_raw_counts)
row.names(GSE71046_raw_counts) <- GSE71046_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE71046 <- GSE71046_raw_counts$Length
names(gene_length_GSE71046) <- GSE71046_raw_counts$Geneid
##################################################
GSE71046_raw_counts <- GSE71046_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(GSE71046_raw_counts) <- c("GSM1826055",
                        "GSM1826056",
                        "GSM1826071",
                        "GSM1826072")
GSE71046_grow_CG <- c("control","control",
                       "drought","drought")
GSE71046_grow_CG <- factor(GSE71046_grow_CG, levels = c("control","drought"))
GSE71046_grow_CG <- relevel(GSE71046_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
GSE71046_sample_info <- data.frame(condition=GSE71046_grow_CG,
                          row.names=names(GSE71046_raw_counts))

dds_GSE71046 <- DESeqDataSetFromMatrix(countData = GSE71046_raw_counts,
                                   colData = GSE71046_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
#assay(DESeq.ds) %>% head

##################################
### filter out count datas and normalization
################################## 
dds_GSE71046 <- dds_GSE71046[rowSums(counts(dds_GSE71046))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE71046)) # number of genes counted
head(assay(dds_GSE71046))

######normalization
dds_GSE71046 <- estimateSizeFactors(dds_GSE71046)
sizeFactors(dds_GSE71046)
sfn_dds_GSE71046 <- counts(dds_GSE71046, normalized = TRUE)

########transformation
ln_dds_GSE71046 <- log2(sfn_dds_GSE71046+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE71046, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE71046 , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE71046 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE71046 <- vst(dds_GSE71046, blind = TRUE)
vstc_dds_GSE71046 <- assay(vst_dds_GSE71046)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE71046,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE71046,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE71046 <- as.dist(1- cor(vstc_dds_GSE71046, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE71046),
      labels = colnames(dis_vstc_dds_GSE71046),
      xlab = "",
      main = "Hierarchical Clustering of GSE71046")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE71046))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE71046)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE71046)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE71046")
print (P)

############## DEG analysis ############################
dds_GSE71046 <- DESeq(dds_GSE71046)
dds_GSE71046 <- estimateSizeFactors(dds_GSE71046)
dds_GSE71046 <- estimateDispersions(dds_GSE71046)
dds_GSE71046 <- nbinomWaldTest(dds_GSE71046)

dds_GSE71046_results <- results(dds_GSE71046, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE71046_results)

dn_GSE71046 <- rownames(subset(dds_GSE71046_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE71046)
