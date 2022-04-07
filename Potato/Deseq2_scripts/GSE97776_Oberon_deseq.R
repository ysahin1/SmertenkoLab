library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
GSE97776_oberon_raw_counts <- read.table("../counts/GSE97776_oberon_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE97776_oberon_raw_counts)
row.names(GSE97776_oberon_raw_counts) <- GSE97776_oberon_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE97776_oberon <- GSE97776_oberon_raw_counts$Length
names(gene_length_GSE97776_oberon) <- GSE97776_oberon_raw_counts$Geneid
##################################################
GSE97776_oberon_raw_counts <- GSE97776_oberon_raw_counts[,-c(1:6)] #c(1:7),10) exclude all information that do not include counts
names(GSE97776_oberon_raw_counts) <- c("GSM2577280",
                        "GSM2577281",
                        "GSM2577282",
                        "GSM2577283",
                        "GSM2577284_1",
                        "GSM2577284_2",
                        "GSM2577285_1",
                        "GSM2577285_2",
                        "GSM2577286_1",
                        "GSM2577286_2",
                        "GSM2577287",
                        "GSM2577288_1",
                        "GSM2577288_2")
GSE97776_oberon_grow_CG <- c("control","control","control",
                          "drought_1","drought_1","drought_1","drought_1","drought_1",
                          "drought_2","drought_2","drought_2","drought_2","drought_2")
#GSE97776_oberon_grow_CG <- factor(GSE97776_oberon_grow_CG, levels = c("control","drought"))
#GSE97776_oberon_grow_CG <- relevel(GSE97776_oberon_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
head(GSE97776_oberon_raw_counts)
GSE97776_oberon_sample_info <- data.frame(condition=GSE97776_oberon_grow_CG,
                          row.names=names(GSE97776_oberon_raw_counts))

dds_GSE97776_oberon <- DESeqDataSetFromMatrix(countData = GSE97776_oberon_raw_counts,
                                   colData = GSE97776_oberon_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
dim(assay(dds_GSE97776_oberon))

##################################
### filter out count datas and normalization
################################## 
dds_GSE97776_oberon <- dds_GSE97776_oberon[rowSums(counts(dds_GSE97776_oberon))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE97776_oberon)) # number of genes counted
head(assay(dds_GSE97776_oberon))

######normalization
dds_GSE97776_oberon <- estimateSizeFactors(dds_GSE97776_oberon)
sizeFactors(dds_GSE97776_oberon)
sfn_dds_GSE97776_oberon <- counts(dds_GSE97776_oberon, normalized = TRUE)

########transformation
ln_dds_GSE97776_oberon <- log2(sfn_dds_GSE97776_oberon+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE97776_oberon, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE97776_oberon , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE97776_oberon [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE97776_oberon <- vst(dds_GSE97776_oberon, blind = TRUE)
vstc_dds_GSE97776_oberon <- assay(vst_dds_GSE97776_oberon)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE97776_oberon,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE97776_oberon,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical cluste  ring##########

dis_vstc_dds_GSE97776_oberon <- as.dist(1- cor(vstc_dds_GSE97776_oberon, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE97776_oberon),
      labels = colnames(dis_vstc_dds_GSE97776_oberon),
      xlab = "",
      main = "Hierarchical Clustering of GSE97776_oberon")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE97776_oberon))
plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_GSE97776_oberon)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE97776_oberon)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE97776_oberon") + geom_text(aes(label=names(GSE97776_oberon_raw_counts)))
print (P)

############## DEG analysis ############################
dds_GSE97776_oberon <- DESeq(dds_GSE97776_oberon)
dds_GSE97776_oberon <- estimateSizeFactors(dds_GSE97776_oberon)
dds_GSE97776_oberon <- estimateDispersions(dds_GSE97776_oberon)
dds_GSE97776_oberon <- nbinomWaldTest(dds_GSE97776_oberon)

dds_GSE97776_oberon_results <- results(dds_GSE97776_oberon, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE97776_oberon_results)

dn_GSE97776_oberon <- rownames(subset(dds_GSE97776_oberon_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE97776_oberon)
dn_GSE97776_oberon[c("PGSC0003DMG400018122","PGSC0003DMG400027073")]
GSE97776_oberon_raw_counts[c("PGSC0003DMG400018122","PGSC0003DMG400027073"),]
