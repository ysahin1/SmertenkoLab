library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
GSE80699_raw_counts <- read.table("GSE80699_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE80699_raw_counts)
row.names(GSE80699_raw_counts) <- GSE80699_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE80699 <- GSE80699_raw_counts$Length
names(gene_length_GSE80699) <- GSE80699_raw_counts$Geneid
##################################################
GSE80699_raw_counts <- GSE80699_raw_counts[,-c(1:6)] #c(1:7),10) exclude all information that do not include counts
names(GSE80699_raw_counts) <- c("GSM2133750",
                        "GSM2133751",
                        "GSM2133752",
                        "GSM2133753",
                        "GSM2133754",
                        "GSM2133755")
GSE80699_grow_CG <- c("control","control","control","drought","drought","drought")
#GSE80699_grow_CG <- factor(GSE80699_grow_CG, levels = c("control","drought"))
#GSE80699_grow_CG <- relevel(GSE80699_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
head(GSE80699_raw_counts)
GSE80699_sample_info <- data.frame(condition=GSE80699_grow_CG,
                          row.names=names(GSE80699_raw_counts))

dds_GSE80699 <- DESeqDataSetFromMatrix(countData = GSE80699_raw_counts,
                                   colData = GSE80699_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
dim(assay(dds_GSE80699))

##################################
### filter out count datas and normalization
################################## 
dds_GSE80699 <- dds_GSE80699[rowSums(counts(dds_GSE80699))> 0,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE80699)) # number of genes counted
head(assay(dds_GSE80699))

######normalization
dds_GSE80699 <- estimateSizeFactors(dds_GSE80699)
sizeFactors(dds_GSE80699)
sfn_dds_GSE80699 <- counts(dds_GSE80699, normalized = TRUE)

########transformation
ln_dds_GSE80699 <- log2(sfn_dds_GSE80699+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE80699, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE80699 , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE80699 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE80699 <- vst(dds_GSE80699, blind = TRUE)
vstc_dds_GSE80699 <- assay(vst_dds_GSE80699)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE80699,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE80699,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE80699 <- as.dist(1- cor(vstc_dds_GSE80699, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE80699),
      labels = colnames(dis_vstc_dds_GSE80699),
      xlab = "",
      main = "Hierarchical Clustering of GSE80699")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE80699))
plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE80699)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE80699)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE80699") + geom_text(aes(label=names(GSE80699_raw_counts)))
print (P)

############## DEG analysis ############################
dds_GSE80699 <- DESeq(dds_GSE80699)
dds_GSE80699 <- estimateSizeFactors(dds_GSE80699)
dds_GSE80699 <- estimateDispersions(dds_GSE80699)
dds_GSE80699 <- nbinomWaldTest(dds_GSE80699)

dds_GSE80699_results <- results(dds_GSE80699, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE80699_results)

dn_GSE80699 <- rownames(subset(dds_GSE80699_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE80699)
dn_GSE80699[c("Sobic.004G200200","Sobic.006G119900")]
