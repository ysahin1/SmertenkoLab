library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
GSE137780_W22_raw_counts <- read.table("../GSE137780_W22_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE137780_W22_raw_counts)
row.names(GSE137780_W22_raw_counts) <- GSE137780_W22_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE137780_W22 <- GSE137780_W22_raw_counts$Length
names(gene_length_GSE137780_W22) <- GSE137780_W22_raw_counts$Geneid
##################################################
GSE137780_W22_raw_counts <- GSE137780_W22_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(GSE137780_W22_raw_counts) <- c("SRR10153108",
                        "SRR10153109",
                        "SRR10153120",
                        "SRR10153121")
GSE137780_W22_grow_CG <- c("control","control","drought","drought")
#GSE137780_W22_grow_CG <- factor(GSE137780_W22_grow_CG, levels = c("0","1","2"))
#GSE137780_W22_grow_CG <- relevel(GSE137780_W22_grow_CG, ref = "0")
####################################
### generate times series deseq object
####################################
GSE137780_W22_sample_info <- data.frame(condition=GSE137780_W22_grow_CG,
                          row.names=names(GSE137780_W22_raw_counts))

dds_GSE137780_W22 <- DESeqDataSetFromMatrix(countData = GSE137780_W22_raw_counts,
                                   colData = GSE137780_W22_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
#assay(DESeq.ds) %>% head

##################################
### filter out count datas and normalization
################################## 
dds_GSE137780_W22 <- dds_GSE137780_W22[rowSums(counts(dds_GSE137780_W22))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE137780_W22)) # number of genes counted
head(assay(dds_GSE137780_W22))

######normalization
dds_GSE137780_W22 <- estimateSizeFactors(dds_GSE137780_W22)
sizeFactors(dds_GSE137780_W22)
sfn_dds_GSE137780_W22 <- counts(dds_GSE137780_W22, normalized = TRUE)

########transformation
ln_dds_GSE137780_W22 <- log2(sfn_dds_GSE137780_W22+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE137780_W22, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE137780_W22 , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE137780_W22 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE137780_W22 <- vst(dds_GSE137780_W22, blind = TRUE)
vstc_dds_GSE137780_W22 <- assay(vst_dds_GSE137780_W22)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE137780_W22,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE137780_W22,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE137780_W22 <- as.dist(1- cor(vstc_dds_GSE137780_W22, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE137780_W22),
      labels = colnames(dis_vstc_dds_GSE137780_W22),
      xlab = "",
      main = "Hierarchical Clustering of GSE137780_W22")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE137780_W22))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE137780_W22)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE137780_W22)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE137780_W22")
print (P)

############## DEG analysis ############################
dds_GSE137780_W22 <- DESeq(dds_GSE137780_W22)
dds_GSE137780_W22 <- estimateSizeFactors(dds_GSE137780_W22)
dds_GSE137780_W22 <- estimateDispersions(dds_GSE137780_W22)
dds_GSE137780_W22 <- nbinomWaldTest(dds_GSE137780_W22)

dds_GSE137780_W22_results <- results(dds_GSE137780_W22, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE137780_W22_results)

dn_GSE137780_W22 <- rownames(subset(dds_GSE137780_W22_results, 
                             padj < 0.05))

length(dn_GSE137780_W22)

dn_GSE137780_W22 <- rownames(subset(dds_GSE137780_W22_results, 
                                    padj < 0.05 & abs(log2FoldChange) > 0))
dn_GSE137780_W22[c("Zm00001d017072", 
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

GSE137780_W22_raw_counts[c("Zm00001d017072", 
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
