library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
GSE97776_Gwiazda_raw_counts <- read.table("../counts/GSE97776_Gwiazda_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE97776_Gwiazda_raw_counts)
row.names(GSE97776_Gwiazda_raw_counts) <- GSE97776_Gwiazda_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE97776_Gwiazda <- GSE97776_Gwiazda_raw_counts$Length
names(gene_length_GSE97776_Gwiazda) <- GSE97776_Gwiazda_raw_counts$Geneid
##################################################
GSE97776_Gwiazda_raw_counts <- GSE97776_Gwiazda_raw_counts[,-c(1:6)] #c(1:7),10) exclude all information that do not include counts
names(GSE97776_Gwiazda_raw_counts) <- c("GSM2577271",
                        "GSM2577272",
                        "GSM2577273",
                        "GSM2577274",
                        "GSM2577275",
                        "GSM2577276",
                        "GSM2577277",
                        "GSM2577278",
                        "GSM2577279")
GSE97776_Gwiazda_grow_CG <- c("control","control","control",
                          "drought_1","drought_1","drought_1",
                          "drought_2","drought_2","drought_2")
#GSE97776_Gwiazda_grow_CG <- factor(GSE97776_Gwiazda_grow_CG, levels = c("control","drought"))
#GSE97776_Gwiazda_grow_CG <- relevel(GSE97776_Gwiazda_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
head(GSE97776_Gwiazda_raw_counts)
GSE97776_Gwiazda_sample_info <- data.frame(condition=GSE97776_Gwiazda_grow_CG,
                          row.names=names(GSE97776_Gwiazda_raw_counts))

dds_GSE97776_Gwiazda <- DESeqDataSetFromMatrix(countData = GSE97776_Gwiazda_raw_counts,
                                   colData = GSE97776_Gwiazda_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
dim(assay(dds_GSE97776_Gwiazda))

##################################
### filter out count datas and normalization
################################## 
dds_GSE97776_Gwiazda <- dds_GSE97776_Gwiazda[rowSums(counts(dds_GSE97776_Gwiazda))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE97776_Gwiazda)) # number of genes counted
head(assay(dds_GSE97776_Gwiazda))

######normalization
dds_GSE97776_Gwiazda <- estimateSizeFactors(dds_GSE97776_Gwiazda)
sizeFactors(dds_GSE97776_Gwiazda)
sfn_dds_GSE97776_Gwiazda <- counts(dds_GSE97776_Gwiazda, normalized = TRUE)

########transformation
ln_dds_GSE97776_Gwiazda <- log2(sfn_dds_GSE97776_Gwiazda+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE97776_Gwiazda, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE97776_Gwiazda , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE97776_Gwiazda [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE97776_Gwiazda <- vst(dds_GSE97776_Gwiazda, blind = TRUE)
vstc_dds_GSE97776_Gwiazda <- assay(vst_dds_GSE97776_Gwiazda)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE97776_Gwiazda,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE97776_Gwiazda,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE97776_Gwiazda <- as.dist(1- cor(vstc_dds_GSE97776_Gwiazda, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE97776_Gwiazda),
      labels = colnames(dis_vstc_dds_GSE97776_Gwiazda),
      xlab = "",
      main = "Hierarchical Clustering of GSE97776_Gwiazda")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE97776_Gwiazda))
plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE97776_Gwiazda)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE97776_Gwiazda)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE97776_Gwiazda") + geom_text(aes(label=names(GSE97776_Gwiazda_raw_counts)))
print (P)

############## DEG analysis ############################
dds_GSE97776_Gwiazda <- DESeq(dds_GSE97776_Gwiazda)
dds_GSE97776_Gwiazda <- estimateSizeFactors(dds_GSE97776_Gwiazda)
dds_GSE97776_Gwiazda <- estimateDispersions(dds_GSE97776_Gwiazda)
dds_GSE97776_Gwiazda <- nbinomWaldTest(dds_GSE97776_Gwiazda)

dds_GSE97776_Gwiazda_results <- results(dds_GSE97776_Gwiazda, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE97776_Gwiazda_results)

dn_GSE97776_Gwiazda <- rownames(subset(dds_GSE97776_Gwiazda_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE97776_Gwiazda)
dn_GSE97776_Gwiazda[c("PGSC0003DMG400018122","PGSC0003DMG400027073")]
GSE97776_Gwiazda_raw_counts[c("PGSC0003DMG400018122","PGSC0003DMG400027073"),]
