library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
GSE48507_raw_counts <- read.table("../GSE48507_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE48507_raw_counts)
row.names(GSE48507_raw_counts) <- GSE48507_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE48507 <- GSE48507_raw_counts$Length
names(gene_length_GSE48507) <- GSE48507_raw_counts$Geneid
##################################################
GSE48507_raw_counts <- GSE48507_raw_counts[,-c(1:6,12)] #exclude all information that do not include counts
names(GSE48507_raw_counts) <- c("GSM1180065",
                        "GSM1180066",
                        "GSM1180067",
                        "GSM1180068",
                        "GSM1180069",
                        "GSM1180070")
GSE48507_grow_CG <- c("0","0","1","1","2","2")
GSE48507_grow_CG <- factor(GSE48507_grow_CG, levels = c("0","1","2"))
GSE48507_grow_CG <- relevel(GSE48507_grow_CG, ref = "0")
####################################
### generate times series deseq object
####################################
GSE48507_sample_info <- data.frame(condition=GSE48507_grow_CG,
                          row.names=names(GSE48507_raw_counts))

dds_GSE48507 <- DESeqDataSetFromMatrix(countData = GSE48507_raw_counts,
                                   colData = GSE48507_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
#assay(DESeq.ds) %>% head

##################################
### filter out count datas and normalization
################################## 
dds_GSE48507 <- dds_GSE48507[rowSums(counts(dds_GSE48507))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE48507)) # number of genes counted
head(assay(dds_GSE48507))

######normalization
dds_GSE48507 <- estimateSizeFactors(dds_GSE48507)
sizeFactors(dds_GSE48507)
sfn_dds_GSE48507 <- counts(dds_GSE48507, normalized = TRUE)

########transformation
ln_dds_GSE48507 <- log2(sfn_dds_GSE48507+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE48507, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE48507 , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE48507 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE48507 <- vst(dds_GSE48507, blind = TRUE)
vstc_dds_GSE48507 <- assay(vst_dds_GSE48507)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE48507,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE48507,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE48507 <- as.dist(1- cor(vstc_dds_GSE48507, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE48507),
      labels = colnames(dis_vstc_dds_GSE48507),
      xlab = "",
      main = "Hierarchical Clustering of GSE48507")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE48507))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE48507)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE48507)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE48507")
print (P)

############## DEG analysis ############################
dds_GSE48507 <- DESeq(dds_GSE48507)
dds_GSE48507 <- estimateSizeFactors(dds_GSE48507)
dds_GSE48507 <- estimateDispersions(dds_GSE48507)
dds_GSE48507 <- nbinomWaldTest(dds_GSE48507)

dds_GSE48507_results <- results(dds_GSE48507, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE48507_results)

dn_GSE48507 <- rownames(subset(dds_GSE48507_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE48507)
dn_GSE48507[c("Zm00001d017072", "Zm00001d003300")]
dn_GSE48507[c("Zm00001d017072", 
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

GSE48507_raw_counts[c("Zm00001d017072", 
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
