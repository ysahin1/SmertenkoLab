library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
SRP132035_raw_counts <- read.table("./SRP132035_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(SRP132035_raw_counts)
row.names(SRP132035_raw_counts) <- SRP132035_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_SRP132035 <- SRP132035_raw_counts$Length
names(gene_length_SRP132035) <- SRP132035_raw_counts$Geneid
##################################################
SRP132035_raw_counts <- SRP132035_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(SRP132035_raw_counts) <- c("SRR6665368",
                        "SRR6665369",
                        "SRR6665370",
                        "SRR6665371",
                        "SRR6665372",
                        "SRR6665373",
                        "SRR6665374",
                        "SRR6665375",
                        "SRR6665376",
                        "SRR6665377",
                        "SRR6665378",
                        "SRR6665379")
SRP132035_grow_CG <- c("control","control",
                       "drought","drought","drought","drought","drought","drought",
                       "control","control","control","control")
SRP132035_grow_CG <- factor(SRP132035_grow_CG, levels = c("control","drought"))
SRP132035_grow_CG <- relevel(SRP132035_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
SRP132035_sample_info <- data.frame(condition=SRP132035_grow_CG,
                          row.names=names(SRP132035_raw_counts))

dds_SRP132035 <- DESeqDataSetFromMatrix(countData = SRP132035_raw_counts,
                                   colData = SRP132035_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
#assay(DESeq.ds) %>% head

##################################
### filter out count datas and normalization
################################## 
dds_SRP132035 <- dds_SRP132035[rowSums(counts(dds_SRP132035))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_SRP132035)) # number of genes counted
head(assay(dds_SRP132035))

######normalization
dds_SRP132035 <- estimateSizeFactors(dds_SRP132035)
sizeFactors(dds_SRP132035)
sfn_dds_SRP132035 <- counts(dds_SRP132035, normalized = TRUE)

########transformation
ln_dds_SRP132035 <- log2(sfn_dds_SRP132035+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_SRP132035, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_SRP132035 , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_SRP132035 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_SRP132035 <- vst(dds_SRP132035, blind = TRUE)
vstc_dds_SRP132035 <- assay(vst_dds_SRP132035)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_SRP132035,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_SRP132035,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_SRP132035 <- as.dist(1- cor(vstc_dds_SRP132035, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_SRP132035),
      labels = colnames(dis_vstc_dds_SRP132035),
      xlab = "",
      main = "Hierarchical Clustering of SRP132035")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_SRP132035))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_SRP132035)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_SRP132035)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of SRP132035")
print (P)

############## DEG analysis ############################
dds_SRP132035 <- DESeq(dds_SRP132035)
dds_SRP132035 <- estimateSizeFactors(dds_SRP132035)
dds_SRP132035 <- estimateDispersions(dds_SRP132035)
dds_SRP132035 <- nbinomWaldTest(dds_SRP132035)

dds_SRP132035_results <- results(dds_SRP132035, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_SRP132035_results)

dn_SRP132035 <- rownames(subset(dds_SRP132035_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_SRP132035)
dn_SRP132035[c("Zm00001d017072", "Zm00001d003300")]
SRP132035_raw_counts[c("Zm00001d017072","Zm00001d003300"),]
