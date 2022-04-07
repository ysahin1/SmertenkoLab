library(magrittr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(DESeq2)


BiocManager::install(c("HTSFilter"))




###################################################
## load file
###################################################
SRP144573_raw_counts <- read.table("./SRP144573_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(SRP144573_raw_counts)
row.names(SRP144573_raw_counts) <- SRP144573_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_SRP144573 <- SRP144573_raw_counts$Length
names(gene_length_SRP144573) <- SRP144573_raw_counts$Geneid
##################################################
SRP144573_raw_counts <- SRP144573_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(SRP144573_raw_counts) <- c("SRR7110724",
                        "SRR7110725")
SRP144573_grow_CG <- c("control","drought")
#SRP144573_grow_CG <- factor(SRP144573_grow_CG, levels = c("control","drought"))
#SRP144573_grow_CG <- relevel(SRP144573_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
SRP144573_sample_info <- data.frame(condition=SRP144573_grow_CG,
                          row.names=names(SRP144573_raw_counts))

dds_SRP144573 <- DESeqDataSetFromMatrix(countData = SRP144573_raw_counts,
                                   colData = SRP144573_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
#assay(DESeq.ds) %>% head

##################################
### filter out count datas and normalization
################################## 
dds_SRP144573 <- dds_SRP144573[rowSums(counts(dds_SRP144573))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_SRP144573)) # number of genes counted
head(assay(dds_SRP144573))

######normalization
dds_SRP144573 <- estimateSizeFactors(dds_SRP144573)
sizeFactors(dds_SRP144573)
sfn_dds_SRP144573 <- counts(dds_SRP144573, normalized = TRUE)

########transformation
ln_dds_SRP144573 <- log2(sfn_dds_SRP144573+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_SRP144573, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_SRP144573 , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_SRP144573 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_SRP144573 <- vst(dds_SRP144573, blind = TRUE)
vstc_dds_SRP144573 <- assay(vst_dds_SRP144573)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_SRP144573,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_SRP144573,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_SRP144573 <- as.dist(1- cor(vstc_dds_SRP144573, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_SRP144573),
      labels = colnames(dis_vstc_dds_SRP144573),
      xlab = "",
      main = "Hierarchical Clustering of SRP144573")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_SRP144573))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_SRP144573)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_SRP144573)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of SRP144573")
print (P)

############## DEG analysis ############################
dds_SRP144573 <- DESeq(dds_SRP144573)
dds_SRP144573 <- estimateSizeFactors(dds_SRP144573)
dds_SRP144573 <- estimateDispersions(dds_SRP144573)
dds_SRP144573 <- nbinomWaldTest(dds_SRP144573)

dds_SRP144573_results <- results(dds_SRP144573, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_SRP144573_results)

dn_SRP144573 <- rownames(subset(dds_SRP144573_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_SRP144573)


# 6 vs 9 degs
#resultsNames(DESeq.ds)


      