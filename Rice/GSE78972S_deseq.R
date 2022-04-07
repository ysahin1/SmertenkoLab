library(ggplot2)
library(DESeq2)

###################################################
## load file
###################################################
GSE78972s_raw_counts <- read.table("raw_counts/GSE78972_short_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE78972s_raw_counts)
row.names(GSE78972s_raw_counts) <- GSE78972s_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE78972s <- GSE78972s_raw_counts$Length
names(gene_length_GSE78972s) <- GSE78972s_raw_counts$Geneid
##################################################
GSE78972s_raw_counts <- GSE78972s_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(GSE78972s_raw_counts) <- c("GSM2082861",
                        "GSM2082862",
                        "GSM2082863",
                        "GSM2082864",
                        "GSM2082867",
                        "GSM2082868")
GSE78972s_grow_CG <- c("control","control","drought","drought","control","drought")

####################################
### generate times series deseq object
####################################
GSE78972s_sample_info <- data.frame(condition=GSE78972s_grow_CG,
                          row.names=names(GSE78972s_raw_counts))

dds_GSE78972s <- DESeqDataSetFromMatrix(countData = GSE78972s_raw_counts,
                                   colData = GSE78972s_sample_info,
                                   design = ~ condition)


##################################
### filter out count datas and normalization
################################## 
dds_GSE78972s <- dds_GSE78972s[rowSums(counts(dds_GSE78972s))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE78972s)) # number of genes counted
head(assay(dds_GSE78972s))

######normalization
dds_GSE78972s <- estimateSizeFactors(dds_GSE78972s)
sizeFactors(dds_GSE78972s)
sfn_dds_GSE78972s <- counts(dds_GSE78972s, normalized = TRUE)

########transformation
ln_dds_GSE78972s <- log2(sfn_dds_GSE78972s+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE78972s, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE78972s , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE78972s [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE78972s <- vst(dds_GSE78972s, blind = TRUE)
vstc_dds_GSE78972s <- assay(vst_dds_GSE78972s)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE78972s,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE78972s,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE78972s <- as.dist(1- cor(vstc_dds_GSE78972s, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE78972s),
      labels = colnames(dis_vstc_dds_GSE78972s),
      xlab = "",
      main = "Hierarchical Clustering of GSE78972s")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE78972s))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE78972s)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE78972s)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE78972s")
print (P)

############## DEG analysis ############################
dds_GSE78972s <- DESeq(dds_GSE78972s)
dds_GSE78972s <- estimateSizeFactors(dds_GSE78972s)
dds_GSE78972s <- estimateDispersions(dds_GSE78972s)
dds_GSE78972s <- nbinomWaldTest(dds_GSE78972s)

dds_GSE78972s_results <- results(dds_GSE78972s, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE78972s_results)

dn_GSE78972s <- rownames(subset(dds_GSE78972s_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE78972s)
