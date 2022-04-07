library(ggplot2)
library(DESeq2)

###################################################
## load file
###################################################
GSE132113_raw_counts <- read.table("../GSE132113_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE132113_raw_counts)
row.names(GSE132113_raw_counts) <- GSE132113_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE132113 <- GSE132113_raw_counts$Length
names(gene_length_GSE132113) <- GSE132113_raw_counts$Geneid
##################################################
GSE132113_raw_counts <- GSE132113_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(GSE132113_raw_counts) <- c("GSM3844504",
                        "GSM3844505",
                        "GSM3844511",
                        "GSM3844512")
GSE132113_grow_CG <- c("control","control","drought","drought")
#GSE132113_grow_CG <- factor(GSE132113_grow_CG, levels = c("control","drought"))
#GSE132113_grow_CG <- relevel(GSE132113_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
GSE132113_sample_info <- data.frame(condition=GSE132113_grow_CG,
                          row.names=names(GSE132113_raw_counts))

dds_GSE132113 <- DESeqDataSetFromMatrix(countData = GSE132113_raw_counts,
                                   colData = GSE132113_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
#assay(DESeq.ds) %>% head

##################################
### filter out count datas and normalization
################################## 
dds_GSE132113 <- dds_GSE132113[rowSums(counts(dds_GSE132113))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE132113)) # number of genes counted
head(assay(dds_GSE132113))

######normalization
dds_GSE132113 <- estimateSizeFactors(dds_GSE132113)
sizeFactors(dds_GSE132113)
sfn_dds_GSE132113 <- counts(dds_GSE132113, normalized = TRUE)

########transformation
ln_dds_GSE132113 <- log2(sfn_dds_GSE132113+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE132113, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE132113 , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE132113 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE132113 <- vst(dds_GSE132113, blind = TRUE)
vstc_dds_GSE132113 <- assay(vst_dds_GSE132113)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE132113,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE132113,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE132113 <- as.dist(1- cor(vstc_dds_GSE132113, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE132113),
      labels = colnames(dis_vstc_dds_GSE132113),
      xlab = "",
      main = "Hierarchical Clustering of GSE132113")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE132113))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE132113)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE132113)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE132113")
print (P)

############## DEG analysis ############################
dds_GSE132113 <- DESeq(dds_GSE132113)
dds_GSE132113 <- estimateSizeFactors(dds_GSE132113)
dds_GSE132113 <- estimateDispersions(dds_GSE132113)
dds_GSE132113 <- nbinomWaldTest(dds_GSE132113)

dds_GSE132113_results <- results(dds_GSE132113, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE132113_results)

dn_GSE132113 <- rownames(subset(dds_GSE132113_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE132113)
dn_GSE132113[c("Zm00001d017072", "Zm00001d003300")]
GSE132113_raw_counts[c("Zm00001d017072", "Zm00001d003300"),]
dn_GSE132113[c("Zm00001d017072", 
              "Zm00001d003300",
              "Zm00001d030960", 
              "Zm00001d014786", 
              "Zm00001d002084", 
              "Zm00001d041395",
              "Zm00001d003040", 
              "Zm00001d015566",
              "Zm00001d051075",
              "Zm00001d017078",
              "Zm00001d038667",
              "Zm00001d010716",
              "Zm00001d036201",
              "Zm00001d038655", 
              "Zm00001d017072"
)]

GSE132113_raw_counts[c("Zm00001d017072", 
                      "Zm00001d003300",
                      "Zm00001d030960", 
                      "Zm00001d014786", 
                      "Zm00001d002084", 
                      "Zm00001d041395",
                      "Zm00001d003040", 
                      "Zm00001d015566",
                      "Zm00001d051075",
                      "Zm00001d017078",
                      "Zm00001d038667",
                      "Zm00001d010716",
                      "Zm00001d036201",
                      "Zm00001d038655", 
                      "Zm00001d017072"
),]
