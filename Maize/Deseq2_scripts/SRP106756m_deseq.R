library(magrittr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(DESeq2)


###################################################
## load file
###################################################
SRP106756_raw_counts <- read.table("./SRP106756_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(SRP106756_raw_counts)
row.names(SRP106756_raw_counts) <- SRP106756_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_SRP106756 <- SRP106756_raw_counts$Length
names(gene_length_SRP106756) <- SRP106756_raw_counts$Geneid
##################################################
SRP106756_raw_counts <- SRP106756_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(SRP106756_raw_counts) <- c("SRR5520482",
                        "SRR5520483",
                        "SRR5520484",
                        "SRR5520485",
                        "SRR5520486",
                        "SRR5520487",
                        "SRR5520488",
                        "SRR5520489")
SRP106756_grow_CG <- c("0","1","2","3","0","1","2","3")
SRP106756_grow_CG <- factor(SRP106756_grow_CG, levels = c("0","1","2","3"))
SRP106756_grow_CG <- relevel(SRP106756_grow_CG, ref = "0")
####################################
### generate times series deseq object
####################################
SRP106756_sample_info <- data.frame(condition=SRP106756_grow_CG,
                          row.names=names(SRP106756_raw_counts))

dds_SRP106756 <- DESeqDataSetFromMatrix(countData = SRP106756_raw_counts,
                                   colData = SRP106756_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
#assay(DESeq.ds) %>% head

##################################
### filter out count datas and normalization
################################## 
dds_SRP106756 <- dds_SRP106756[rowSums(counts(dds_SRP106756))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_SRP106756)) # number of genes counted
head(assay(dds_SRP106756))

######normalization
dds_SRP106756 <- estimateSizeFactors(dds_SRP106756)
sizeFactors(dds_SRP106756)
sfn_dds_SRP106756 <- counts(dds_SRP106756, normalized = TRUE)

########transformation
ln_dds_SRP106756 <- log2(sfn_dds_SRP106756+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_SRP106756, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_SRP106756 , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_SRP106756 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_SRP106756 <- vst(dds_SRP106756, blind = TRUE)
vstc_dds_SRP106756 <- assay(vst_dds_SRP106756)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_SRP106756,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_SRP106756,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_SRP106756 <- as.dist(1- cor(vstc_dds_SRP106756, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_SRP106756),
      labels = colnames(dis_vstc_dds_SRP106756),
      xlab = "",
      main = "Hierarchical Clustering of SRP106756")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_SRP106756))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_SRP106756)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_SRP106756)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of SRP106756")
print (P)

############## DEG analysis ############################
dds_SRP106756 <- DESeq(dds_SRP106756)
dds_SRP106756 <- estimateSizeFactors(dds_SRP106756)
dds_SRP106756 <- estimateDispersions(dds_SRP106756)
dds_SRP106756 <- nbinomWaldTest(dds_SRP106756)

dds_SRP106756_results <- results(dds_SRP106756, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_SRP106756_results)

dn_SRP106756 <- rownames(subset(dds_SRP106756_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_SRP106756)
