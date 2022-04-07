library(ggplot2)
library(DESeq2)
###################################################
## load file
###################################################
GSE151277_d1_raw_counts <- read.table("../counts/GSE151277_D_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(GSE151277_d1_raw_counts)
row.names(GSE151277_d1_raw_counts) <- GSE151277_d1_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_GSE151277_d1 <- GSE151277_d1_raw_counts$Length
names(gene_length_GSE151277_d1) <- GSE151277_d1_raw_counts$Geneid
##################################################
GSE151277_d1_raw_counts <- GSE151277_d1_raw_counts[,-c(1:6,13:24)] #c(1:7),10) exclude all information that do not include counts
names(GSE151277_d1_raw_counts) <- c("GSM4570167",
                                   "GSM4570168",
                                   "GSM4570169",
                                   "GSM4570170",
                                   "GSM4570171",
                                   "GSM4570172")
                                   #"GSM4570173",
                                   #"GSM4570174",
                                   #"GSM4570175",
                                   #"GSM4570176",
                                   #"GSM4570177",
                                   #"GSM4570178",
                                   #"GSM4570179",
                                   #"GSM4570180",
                                   #"GSM4570181",
                                   #"GSM4570182",
                                   #"GSM4570183",
                                   #"GSM4570184")
GSE151277_d1_grow_CG <- c("control","control","control",
                         "drought_1","drought_1","drought_1")
                         #"drought_2","drought_2","drought_2",
                         #"drought_3","drought_3","drought_3",
                         #"drought_4","drought_4","drought_4",
                         #"drought_5","drought_5","drought_5")
#GSE151277_d1_grow_CG <- factor(GSE151277_d1_grow_CG, levels = c("control","drought"))
#GSE151277_d1_grow_CG <- relevel(GSE151277_d1_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
head(GSE151277_d1_raw_counts)
GSE151277_d1_sample_info <- data.frame(condition=GSE151277_d1_grow_CG,
                                      row.names=names(GSE151277_d1_raw_counts))

dds_GSE151277_d1 <- DESeqDataSetFromMatrix(countData = GSE151277_d1_raw_counts,
                                          colData = GSE151277_d1_sample_info,
                                          design = ~ condition)

#colData(DESeq.ds)
dim(assay(dds_GSE151277_d1))

##################################
### filter out count datas and normalization
################################## 
dds_GSE151277_d1 <- dds_GSE151277_d1[rowSums(counts(dds_GSE151277_d1))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_GSE151277_d1)) # number of genes counted
head(assay(dds_GSE151277_d1))

######normalization
dds_GSE151277_d1 <- estimateSizeFactors(dds_GSE151277_d1)
sizeFactors(dds_GSE151277_d1)
sfn_dds_GSE151277_d1 <- counts(dds_GSE151277_d1, normalized = TRUE)

########transformation
ln_dds_GSE151277_d1 <- log2(sfn_dds_GSE151277_d1+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_GSE151277_d1, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_GSE151277_d1 , notch = TRUE ,
        main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
        ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_GSE151277_d1 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_GSE151277_d1 <- vst(dds_GSE151277_d1, blind = TRUE)
vstc_dds_GSE151277_d1 <- assay(vst_dds_GSE151277_d1)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_GSE151277_d1,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_GSE151277_d1,
                       ranks = FALSE, # show the data on the original scale
                       plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_GSE151277_d1 <- as.dist(1- cor(vstc_dds_GSE151277_d1, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_GSE151277_d1),
     labels = colnames(dis_vstc_dds_GSE151277_d1),
     xlab = "",
     main = "Hierarchical Clustering of GSE151277_d1")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_GSE151277_d1))
plot(pc$x[ ,1],pc$x[ ,2],
     col = colData(dds_GSE151277_d1)[ ,1],
     main = "Principal Component Analysis")


P <- plotPCA(vst_dds_GSE151277_d1)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of GSE151277_d1") + geom_text(aes(label=names(GSE151277_d1_raw_counts)))
print (P)

############## DEG analysis ############################
dds_GSE151277_d1 <- DESeq(dds_GSE151277_d1)
dds_GSE151277_d1 <- estimateSizeFactors(dds_GSE151277_d1)
dds_GSE151277_d1 <- estimateDispersions(dds_GSE151277_d1)
dds_GSE151277_d1 <- nbinomWaldTest(dds_GSE151277_d1)

dds_GSE151277_d1_results <- results(dds_GSE151277_d1, 
                                   independentFiltering = TRUE, 
                                   alpha = 0.05, lfcThreshold = 0)
summary(dds_GSE151277_d1_results)

dn_GSE151277_d1 <- rownames(subset(dds_GSE151277_d1_results, 
                                  padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_GSE151277_d1)


  