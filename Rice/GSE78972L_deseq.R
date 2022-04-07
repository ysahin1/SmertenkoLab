library(ggplot2)
library(DESeq2)

###################################################
## load file
###################################################
GSE78972l_raw_counts <- read.table("../GSE78972_long_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE78972l_raw_counts)
row.names(GSE78972l_raw_counts) <- GSE78972l_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE78972l <- GSE78972l_raw_counts$Length
names(gene_length_GSE78972l) <- GSE78972l_raw_counts$Geneid
##################################################
GSE78972l_raw_counts <- GSE78972l_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(GSE78972l_raw_counts) <- c("GSM2082857",
                        "GSM2082858",
                        "GSM2082859",
                        "GSM2082860",
                        "GSM2082865",
                        "GSM2082866")
GSE78972l_grow_CG <- c("control","control","drought","drought","control","drought")

####################################
### generate times series deseq object
####################################
GSE78972l_sample_info <- data.frame(condition=GSE78972l_grow_CG,
                          row.names=names(GSE78972l_raw_counts))

dds_GSE78972l <- DESeqDataSetFromMatrix(countData = GSE78972l_raw_counts,
                                   colData = GSE78972l_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
#assay(DESeq.ds) %>% head

##################################
### filter out count datas and normalization
################################## 
dds_GSE78972l <- dds_GSE78972l[rowSums(counts(dds_GSE78972l))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE78972l)) # number of genes counted
head(assay(dds_GSE78972l))

######normalization
dds_GSE78972l <- estimateSizeFactors(dds_GSE78972l)
sizeFactors(dds_GSE78972l)
sfn_dds_GSE78972l <- counts(dds_GSE78972l, normalized = TRUE)

########transformation
ln_dds_GSE78972l <- log2(sfn_dds_GSE78972l+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE78972l, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE78972l , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE78972l [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE78972l <- vst(dds_GSE78972l, blind = TRUE)
vstc_dds_GSE78972l <- assay(vst_dds_GSE78972l)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE78972l,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE78972l,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE78972l <- as.dist(1- cor(vstc_dds_GSE78972l, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE78972l),
      labels = colnames(dis_vstc_dds_GSE78972l),
      xlab = "",
      main = "Hierarchical Clustering of GSE78972l")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE78972l))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE78972l)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE78972l)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE78972l")
print (P)

############## DEG analysis ############################
dds_GSE78972l <- DESeq(dds_GSE78972l)
dds_GSE78972l <- estimateSizeFactors(dds_GSE78972l)
dds_GSE78972l <- estimateDispersions(dds_GSE78972l)
dds_GSE78972l <- nbinomWaldTest(dds_GSE78972l)

dds_GSE78972l_results <- results(dds_GSE78972l, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE78972l_results)

dn_GSE78972l <- rownames(subset(dds_GSE78972l_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE78972l)
dn_GSE78972l["LOC_Os02g38050"] 