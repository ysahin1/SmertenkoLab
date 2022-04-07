library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
GSE97776_tajfun_raw_counts <- read.table("../counts/GSE97776_tajfun_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE97776_tajfun_raw_counts)
row.names(GSE97776_tajfun_raw_counts) <- GSE97776_tajfun_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE97776_tajfun <- GSE97776_tajfun_raw_counts$Length
names(gene_length_GSE97776_tajfun) <- GSE97776_tajfun_raw_counts$Geneid
##################################################
GSE97776_tajfun_raw_counts <- GSE97776_tajfun_raw_counts[,-c(1:6)] #c(1:7),10) exclude all information that do not include counts
names(GSE97776_tajfun_raw_counts) <- c("GSM2577289_1",
                                       "GSM2577289_2",
                                       "GSM2577290_1",
                                       "GSM2577290_2",
                                       "GSM2577291_1",
                                       "GSM2577291_2",
                                       "GSM2577292_1",
                                       "GSM2577292_2",
                                       "GSM2577293_1",
                                       "GSM2577293_2",
                                       "GSM2577294_1",
                                       "GSM2577294_2",
                                       "GSM2577295_1",
                                       "GSM2577295_2",
                                       "GSM2577296",
                                       "GSM2577297_1",
                                       "GSM2577297_2")
GSE97776_tajfun_grow_CG <- c("control","control","control","control","control","control",
                          "drought_1","drought_1","drought_1","drought_1","drought_1","drought_1",
                          "drought_2","drought_2","drought_2","drought_2","drought_2")
#GSE97776_tajfun_grow_CG <- factor(GSE97776_tajfun_grow_CG, levels = c("control","drought"))
#GSE97776_tajfun_grow_CG <- relevel(GSE97776_tajfun_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
head(GSE97776_tajfun_raw_counts)
GSE97776_tajfun_sample_info <- data.frame(condition=GSE97776_tajfun_grow_CG,
                          row.names=names(GSE97776_tajfun_raw_counts))

dds_GSE97776_tajfun <- DESeqDataSetFromMatrix(countData = GSE97776_tajfun_raw_counts,
                                   colData = GSE97776_tajfun_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
dim(assay(dds_GSE97776_tajfun))

##################################
### filter out count datas and normalization
################################## 
dds_GSE97776_tajfun <- dds_GSE97776_tajfun[rowSums(counts(dds_GSE97776_tajfun))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE97776_tajfun)) # number of genes counted
head(assay(dds_GSE97776_tajfun))

######normalization
dds_GSE97776_tajfun <- estimateSizeFactors(dds_GSE97776_tajfun)
sizeFactors(dds_GSE97776_tajfun)
sfn_dds_GSE97776_tajfun <- counts(dds_GSE97776_tajfun, normalized = TRUE)

########transformation
ln_dds_GSE97776_tajfun <- log2(sfn_dds_GSE97776_tajfun+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE97776_tajfun, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE97776_tajfun , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE97776_tajfun [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE97776_tajfun <- vst(dds_GSE97776_tajfun, blind = TRUE)
vstc_dds_GSE97776_tajfun <- assay(vst_dds_GSE97776_tajfun)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE97776_tajfun,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE97776_tajfun,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical cluste  ring##########

dis_vstc_dds_GSE97776_tajfun <- as.dist(1- cor(vstc_dds_GSE97776_tajfun, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE97776_tajfun),
      labels = colnames(dis_vstc_dds_GSE97776_tajfun),
      xlab = "",
      main = "Hierarchical Clustering of GSE97776_tajfun")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE97776_tajfun))
plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE97776_tajfun)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE97776_tajfun)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE97776_tajfun") + geom_text(aes(label=names(GSE97776_tajfun_raw_counts)))
print (P)

############## DEG analysis ############################
dds_GSE97776_tajfun <- DESeq(dds_GSE97776_tajfun)
dds_GSE97776_tajfun <- estimateSizeFactors(dds_GSE97776_tajfun)
dds_GSE97776_tajfun <- estimateDispersions(dds_GSE97776_tajfun)
dds_GSE97776_tajfun <- nbinomWaldTest(dds_GSE97776_tajfun)

dds_GSE97776_tajfun_results <- results(dds_GSE97776_tajfun, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE97776_tajfun_results)

dn_GSE97776_tajfun <- rownames(subset(dds_GSE97776_tajfun_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE97776_tajfun)
dn_GSE97776_tajfun[c("PGSC0003DMG400018122","PGSC0003DMG400027073")]
GSE97776_tajfun_raw_counts[c("PGSC0003DMG400018122","PGSC0003DMG400027073"),]
