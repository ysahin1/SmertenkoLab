library(magrittr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(DESeq2)


###################################################
## load file
###################################################
GSE57950h_raw_counts <- read.table("../H471_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE57950h_raw_counts)
row.names(GSE57950h_raw_counts) <- GSE57950h_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE57950h <- GSE57950h_raw_counts$Length
names(gene_length_GSE57950h) <- GSE57950h_raw_counts$Geneid
##################################################
GSE57950h_raw_counts <- GSE57950h_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(GSE57950h_raw_counts) <- c("GSM1398433",
                        "GSM1398434",
                        "GSM1398445",
                        "GSM1398446")
GSE57950h_grow_CG <- c("control","control","drought","drought")
#GSE57950h_grow_CG <- factor(GSE57950h_grow_CG, levels = c("control","drought"))
#GSE57950h_grow_CG <- relevel(GSE57950h_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
GSE57950h_sample_info <- data.frame(condition=GSE57950h_grow_CG,
                          row.names=names(GSE57950h_raw_counts))

dds_GSE57950h <- DESeqDataSetFromMatrix(countData = GSE57950h_raw_counts,
                                   colData = GSE57950h_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
#assay(DESeq.ds) %>% head

##################################
### filter out count datas and normalization
################################## 
dds_GSE57950h <- dds_GSE57950h[rowSums(counts(dds_GSE57950h))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE57950h)) # number of genes counted
head(assay(dds_GSE57950h))

######normalization
dds_GSE57950h <- estimateSizeFactors(dds_GSE57950h)
sizeFactors(dds_GSE57950h)
sfn_dds_GSE57950h <- counts(dds_GSE57950h, normalized = TRUE)

########transformation
ln_dds_GSE57950h <- log2(sfn_dds_GSE57950h+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE57950h, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE57950h , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE57950h [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE57950h <- vst(dds_GSE57950h, blind = TRUE)
vstc_dds_GSE57950h <- assay(vst_dds_GSE57950h)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE57950h,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE57950h,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE57950h <- as.dist(1- cor(vstc_dds_GSE57950h, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE57950h),
      labels = colnames(dis_vstc_dds_GSE57950h),
      xlab = "",
      main = "Hierarchical Clustering of GSE57950h")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE57950h))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE57950h)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE57950h)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE57950h")
print (P)

############## DEG analysis ############################
dds_GSE57950h <- DESeq(dds_GSE57950h)
dds_GSE57950h <- estimateSizeFactors(dds_GSE57950h)
dds_GSE57950h <- estimateDispersions(dds_GSE57950h)
dds_GSE57950h <- nbinomWaldTest(dds_GSE57950h)

dds_GSE57950h_results <- results(dds_GSE57950h, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE57950h_results)

dn_GSE57950h <- rownames(subset(dds_GSE57950h_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE57950h)
dn_GSE57950h[c("LOC_Os06g05190",
               "LOC_Os09g26870",
               "LOC_Os04g55480",
               "LOC_Os03g49210",
               "LOC_Os05g40810",
               "LOC_Os02g39990",
               "LOC_Os02g09060",
               "LOC_Os03g59650",
               "LOC_Os05g43610",
               "LOC_Os06g05720",
               "LOC_Os06g50170",
               "LOC_Os05g43280",
               "LOC_Os05g46490",
               "LOC_Os08g31930",
               "LOC_Os04g43300",
               "LOC_Os02g38050")]
GSE57950h_raw_counts[c("LOC_Os06g05190",
                       "LOC_Os09g26870",
                       "LOC_Os04g55480",
                       "LOC_Os03g49210",
                       "LOC_Os05g40810",
                       "LOC_Os02g39990",
                       "LOC_Os02g09060",
                       "LOC_Os03g59650",
                       "LOC_Os05g43610",
                       "LOC_Os06g05720",
                       "LOC_Os06g50170",
                       "LOC_Os05g43280",
                       "LOC_Os05g46490",
                       "LOC_Os08g31930",
                       "LOC_Os04g43300",
                       "LOC_Os02g38050"),]
