library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
GSE65022_raw_counts <- read.table("raw_counts/GSE65022_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE65022_raw_counts)
row.names(GSE65022_raw_counts) <- GSE65022_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE65022 <- GSE65022_raw_counts$Length
names(gene_length_GSE65022) <- GSE65022_raw_counts$Geneid
##################################################
GSE65022_raw_counts <- GSE65022_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(GSE65022_raw_counts) <- c("GSM1585980",
                        "GSM1585981",
                        "GSM1585982",
                        "GSM1585983")
GSE65022_grow_CG <- c("drought","drought","control","control")
GSE65022_grow_CG <- factor(GSE65022_grow_CG, levels = c("control","drought"))
GSE65022_grow_CG <- relevel(GSE65022_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
GSE65022_sample_info <- data.frame(condition=GSE65022_grow_CG,
                          row.names=names(GSE65022_raw_counts))

dds_GSE65022 <- DESeqDataSetFromMatrix(countData = GSE65022_raw_counts,
                                   colData = GSE65022_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
#assay(DESeq.ds) %>% head

##################################
### filter out count datas and normalization
################################## 
dds_GSE65022 <- dds_GSE65022[rowSums(counts(dds_GSE65022))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE65022)) # number of genes counted
head(assay(dds_GSE65022))

######normalization
dds_GSE65022 <- estimateSizeFactors(dds_GSE65022)
sizeFactors(dds_GSE65022)
sfn_dds_GSE65022 <- counts(dds_GSE65022, normalized = TRUE)

########transformation
ln_dds_GSE65022 <- log2(sfn_dds_GSE65022+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE65022, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE65022 , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE65022 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE65022 <- vst(dds_GSE65022, blind = TRUE)
vstc_dds_GSE65022 <- assay(vst_dds_GSE65022)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE65022,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE65022,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE65022 <- as.dist(1- cor(vstc_dds_GSE65022, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE65022),
      labels = colnames(dis_vstc_dds_GSE65022),
      xlab = "",
      main = "Hierarchical Clustering of GSE65022")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE65022))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE65022)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE65022)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE65022")
print (P)

############## DEG analysis ############################
dds_GSE65022 <- DESeq(dds_GSE65022)
dds_GSE65022_results <- results(dds_GSE65022, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE65022_results)

dn_GSE65022 <- rownames(subset(dds_GSE65022_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE65022)
