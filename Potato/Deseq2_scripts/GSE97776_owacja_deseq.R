library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
GSE97776_owacja_raw_counts <- read.table("../counts/GSE97776_owacja_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE97776_owacja_raw_counts)
row.names(GSE97776_owacja_raw_counts) <- GSE97776_owacja_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE97776_owacja <- GSE97776_owacja_raw_counts$Length
names(gene_length_GSE97776_owacja) <- GSE97776_owacja_raw_counts$Geneid
##################################################
GSE97776_owacja_raw_counts <- GSE97776_owacja_raw_counts[,-c(1:6)] #c(1:7),10) exclude all information that do not include counts
names(GSE97776_owacja_raw_counts) <- c("GSM2577298_1",
                                       "GSM2577298_2",
                                       "GSM2577299",
                                       "GSM2577300_1",
                                       "GSM2577300_2",
                                       "GSM2577301",
                                       "GSM2577302",
                                       "GSM2577303_1",
                                       "GSM2577303_2",
                                       "GSM2577304_1",
                                       "GSM2577304_2",
                                       "GSM2577305",
                                       "GSM2577306"
                                       )
GSE97776_owacja_grow_CG <- c("control","control","control","control","control",
                          "drought_1","drought_1","drought_1","drought_1",
                          "drought_2","drought_2","drought_2","drought_2")
#GSE97776_owacja_grow_CG <- factor(GSE97776_owacja_grow_CG, levels = c("control","drought"))
#GSE97776_owacja_grow_CG <- relevel(GSE97776_owacja_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
head(GSE97776_owacja_raw_counts)
GSE97776_owacja_sample_info <- data.frame(condition=GSE97776_owacja_grow_CG,
                          row.names=names(GSE97776_owacja_raw_counts))

dds_GSE97776_owacja <- DESeqDataSetFromMatrix(countData = GSE97776_owacja_raw_counts,
                                   colData = GSE97776_owacja_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
dim(assay(dds_GSE97776_owacja))

##################################
### filter out count datas and normalization
################################## 
dds_GSE97776_owacja <- dds_GSE97776_owacja[rowSums(counts(dds_GSE97776_owacja))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE97776_owacja)) # number of genes counted
head(assay(dds_GSE97776_owacja))

######normalization
dds_GSE97776_owacja <- estimateSizeFactors(dds_GSE97776_owacja)
sizeFactors(dds_GSE97776_owacja)
sfn_dds_GSE97776_owacja <- counts(dds_GSE97776_owacja, normalized = TRUE)

########transformation
ln_dds_GSE97776_owacja <- log2(sfn_dds_GSE97776_owacja+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE97776_owacja, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE97776_owacja , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE97776_owacja [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE97776_owacja <- vst(dds_GSE97776_owacja, blind = TRUE)
vstc_dds_GSE97776_owacja <- assay(vst_dds_GSE97776_owacja)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE97776_owacja,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE97776_owacja,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical cluste  ring##########

dis_vstc_dds_GSE97776_owacja <- as.dist(1- cor(vstc_dds_GSE97776_owacja, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE97776_owacja),
      labels = colnames(dis_vstc_dds_GSE97776_owacja),
      xlab = "",
      main = "Hierarchical Clustering of GSE97776_owacja")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE97776_owacja))
plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE97776_owacja)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE97776_owacja)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE97776_owacja") + geom_text(aes(label=names(GSE97776_owacja_raw_counts)))
print (P)

############## DEG analysis ############################
dds_GSE97776_owacja <- DESeq(dds_GSE97776_owacja)
dds_GSE97776_owacja <- estimateSizeFactors(dds_GSE97776_owacja)
dds_GSE97776_owacja <- estimateDispersions(dds_GSE97776_owacja)
dds_GSE97776_owacja <- nbinomWaldTest(dds_GSE97776_owacja)

dds_GSE97776_owacja_results <- results(dds_GSE97776_owacja, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE97776_owacja_results)

dn_GSE97776_owacja <- rownames(subset(dds_GSE97776_owacja_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE97776_owacja)
dn_GSE97776_owacja[c("PGSC0003DMG400018122","PGSC0003DMG400027073")]
GSE97776_owacja_raw_counts[c("PGSC0003DMG400018122","PGSC0003DMG400027073"),]
