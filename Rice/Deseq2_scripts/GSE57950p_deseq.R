library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
GSE57950p_raw_counts <- read.table("raw_counts/P28_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE57950p_raw_counts)
row.names(GSE57950p_raw_counts) <- GSE57950p_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE57950p <- GSE57950p_raw_counts$Length
names(gene_length_GSE57950p) <- GSE57950p_raw_counts$Geneid
##################################################
GSE57950p_raw_counts <- GSE57950p_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(GSE57950p_raw_counts) <- c("GSM1398431",
                        "GSM1398432",
                        "GSM1398443",
                        "GSM1398444")
GSE57950p_grow_CG <- c("control","control","drought","drought")

####################################
### generate times series deseq object
####################################
GSE57950p_sample_info <- data.frame(condition=GSE57950p_grow_CG,
                          row.names=names(GSE57950p_raw_counts))

dds_GSE57950p <- DESeqDataSetFromMatrix(countData = GSE57950p_raw_counts,
                                   colData = GSE57950p_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
#assay(DESeq.ds) %>% head

##################################
### filter out count datas and normalization
################################## 
dds_GSE57950p <- dds_GSE57950p[rowSums(counts(dds_GSE57950p))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE57950p)) # number of genes counted
head(assay(dds_GSE57950p))

######normalization
dds_GSE57950p <- estimateSizeFactors(dds_GSE57950p)
sizeFactors(dds_GSE57950p)
sfn_dds_GSE57950p <- counts(dds_GSE57950p, normalized = TRUE)

########transformation
ln_dds_GSE57950p <- log2(sfn_dds_GSE57950p+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE57950p, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE57950p , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE57950p [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE57950p <- vst(dds_GSE57950p, blind = TRUE)
vstc_dds_GSE57950p <- assay(vst_dds_GSE57950p)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE57950p,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE57950p,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE57950p <- as.dist(1- cor(vstc_dds_GSE57950p, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE57950p),
      labels = colnames(dis_vstc_dds_GSE57950p),
      xlab = "",
      main = "Hierarchical Clustering of GSE57950p")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE57950p))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE57950p)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE57950p)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE57950p")
print (P)

############## DEG analysis ############################
dds_GSE57950p <- DESeq(dds_GSE57950p)
dds_GSE57950p <- estimateSizeFactors(dds_GSE57950p)
dds_GSE57950p <- estimateDispersions(dds_GSE57950p)
dds_GSE57950p <- nbinomWaldTest(dds_GSE57950p)

dds_GSE57950p_results <- results(dds_GSE57950p, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE57950p_results)

dn_GSE57950p <- rownames(subset(dds_GSE57950p_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE57950p)
