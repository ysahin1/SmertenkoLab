library(DESeq2)

GSE132113 <- read.table("../GSE132113_raw_counts_hisat-sorted.txt",
                        header = T)
row.names(GSE132113) <- GSE132113$Geneid

GSE132113 <- GSE132113[,-c(1:6)]
colnames(GSE132113) <- c("Zm.13.C1",
                         "Zm.13.C2",
                         "Zm.13.S1",
                         "Zm.13.S2")
head(GSE132113)
GSE132113_condition <- c("control","control","drought","drought")

GSE132113_sample_info <- data.frame(condition=GSE132113_condition,
                                    row.names=names(GSE132113))

dds_GSE132113 <- DESeqDataSetFromMatrix(countData = GSE132113,
                                        colData = GSE132113_sample_info,
                                        design = ~ condition)
dds_GSE132113 <- dds_GSE132113[rowSums(counts(dds_GSE132113))> 100,]
GSE132113_vst <- vst(dds_GSE132113, blind = F)
GSE132113_vst_counts <- assay(GSE132113_vst)
###########################HHZ########################################
GSE48507 <- read.table("../GSE48507_raw_counts_hisat-sorted.txt",
                         header = T)
row.names(GSE48507) <- GSE48507$Geneid
head(GSE48507)
GSE48507 <- GSE48507[,-c(1:6,12)]
colnames(GSE48507) <- c("Zm.07.C1",
                          "Zm.07.C2",
                          "Zm.07.S1",
                          "Zm.07.S2",
                        "Zm.07.S3",
                        "Zm.07.S4")
head(GSE48507)
GSE48507_condition <- c("control","control","drought","drought","drought","drought")

GSE48507_sample_info <- data.frame(condition=GSE48507_condition,
                                     row.names=names(GSE48507))

dds_GSE48507 <- DESeqDataSetFromMatrix(countData = GSE48507,
                                         colData = GSE48507_sample_info,
                                         design = ~ condition)
dds_GSE48507 <- dds_GSE48507[rowSums(counts(dds_GSE48507))> 100,]
GSE48507_vst <- vst(dds_GSE48507, blind = F)
GSE48507_vst_counts <- assay(GSE48507_vst)

###########################ATH########################################
GSE137780_b104 <- read.table("../GSE137780_B104_raw_counts_hisat-sorted.txt",
                        header = T)
row.names(GSE137780_b104) <- GSE137780_b104$Geneid

GSE137780_b104 <- GSE137780_b104[,-c(1:6)]
colnames(GSE137780_b104) <- c("Zm.b104.80.C1",
                         "Zm.b104.80.C2",
                         "Zm.b104.80.S1",
                         "Zm.b104.80.S2")
head(GSE137780_b104)
GSE137780_b104_condition <- c("control","control","drought","drought")

GSE137780_b104_sample_info <- data.frame(condition=GSE137780_b104_condition,
                                    row.names=names(GSE137780_b104))

dds_GSE137780_b104 <- DESeqDataSetFromMatrix(countData = GSE137780_b104,
                                        colData = GSE137780_b104_sample_info,
                                        design = ~ condition)
dds_GSE137780_b104 <- dds_GSE137780_b104[rowSums(counts(dds_GSE137780_b104))> 100,]
GSE137780_b104_vst <- vst(dds_GSE137780_b104, blind = F)
GSE137780_b104_vst_counts <- assay(GSE137780_b104_vst)

###########################SD########################################
GSE137780_DH4866 <- read.table("../GSE137780_DH4866_raw_counts_hisat-sorted.txt",
                        header = T)
row.names(GSE137780_DH4866) <- GSE137780_DH4866$Geneid

GSE137780_DH4866 <- GSE137780_DH4866[,-c(1:6)]
colnames(GSE137780_DH4866) <- c("Zm.DH4866.80.C1",
                         "Zm.DH4866.80.C2",
                         "Zm.DH4866.80.S1",
                         "Zm.DH4866.80.S2")
head(GSE137780_DH4866)
GSE137780_DH4866_condition <- c("control","control","drought","drought")

GSE137780_DH4866_sample_info <- data.frame(condition=GSE137780_DH4866_condition,
                                    row.names=names(GSE137780_DH4866))

dds_GSE137780_DH4866 <- DESeqDataSetFromMatrix(countData = GSE137780_DH4866,
                                        colData = GSE137780_DH4866_sample_info,
                                        design = ~ condition)
dds_GSE137780_DH4866 <- dds_GSE137780_DH4866[rowSums(counts(dds_GSE137780_DH4866))> 100,]
GSE137780_DH4866_vst <- vst(dds_GSE137780_DH4866, blind = F)
GSE137780_DH4866_vst_counts <- assay(GSE137780_DH4866_vst)

###########################LD########################################
GSE137780_W22 <- read.table("../GSE137780_W22_raw_counts_hisat-sorted.txt",
                        header = T)
row.names(GSE137780_W22) <- GSE137780_W22$Geneid

GSE137780_W22 <- GSE137780_W22[,-c(1:6)]
colnames(GSE137780_W22) <- c("Zm.W22.80.C1",
                         "Zm.W22.80.C2",
                         "Zm.W22.80.S1",
                         "Zm.W22.80.S2")
head(GSE137780_W22)
GSE137780_W22_condition <- c("control","control","drought","drought")

GSE137780_W22_sample_info <- data.frame(condition=GSE137780_W22_condition,
                                    row.names=names(GSE137780_W22))

dds_GSE137780_W22 <- DESeqDataSetFromMatrix(countData = GSE137780_W22,
                                        colData = GSE137780_W22_sample_info,
                                        design = ~ condition)
dds_GSE137780_W22 <- dds_GSE137780_W22[rowSums(counts(dds_GSE137780_W22))> 100,]
GSE137780_W22_vst <- vst(dds_GSE137780_W22, blind = F)
GSE137780_W22_vst_counts <- assay(GSE137780_W22_vst)


# GSE134945 filtering #

GSE132113_rv <- rowVars(GSE132113_vst_counts)
GSE132113_rv <- order(GSE132113_rv, decreasing = TRUE)
GSE132113_mvg <- GSE132113_vst_counts[GSE132113_rv,][1:2000,]

GSE48507_rv <- rowVars(GSE48507_vst_counts)
GSE48507_rv <- order(GSE48507_rv, decreasing = TRUE)
GSE48507_mvg <- GSE48507_vst_counts[GSE48507_rv,][1:2000,]

GSE137780_b104_rv <- rowVars(GSE137780_b104_vst_counts)
GSE137780_b104_rv <- order(GSE137780_b104_rv, decreasing = TRUE)
GSE137780_b104_mvg <- GSE137780_b104_vst_counts[GSE137780_b104_rv,][1:2000,]

GSE137780_DH4866_rv <- rowVars(GSE137780_DH4866_vst_counts)
GSE137780_DH4866_rv <- order(GSE137780_DH4866_rv, decreasing = TRUE)
GSE137780_DH4866_mvg <- GSE137780_DH4866_vst_counts[GSE137780_DH4866_rv,][1:2000,]

GSE137780_W22_rv <- rowVars(GSE137780_W22_vst_counts)
GSE137780_W22_rv <- order(GSE137780_W22_rv, decreasing = TRUE)
GSE137780_W22_mvg <- GSE137780_W22_vst_counts[GSE137780_W22_rv,][1:2000,]


mvg <- as.matrix(cbind(GSE132113_mvg, GSE48507_mvg, GSE137780_b104_mvg, GSE137780_DH4866_mvg, GSE137780_W22_mvg))

dis_mvg <- as.dist(1- cor(mvg, method ="pearson"))
pc <- prcomp (t(dis_mvg))


library(scatterplot3d)
colors <-  c(rep("#999999", 18), 
             rep("#E69F00",20),
             rep("#B2182B", 18))
length(colors)
dim(pc$x[,1:3])
pch <-  c(tm_condition$pch,
          ath_condition$pch,
          zm_condition$pch)
days <- c(tm_condition$days,
          ath_condition$days,
          zm_condition$days)
par(mfrow=c(2,2))
tiff(filename = "3D_PCA.tiff", width = 8, height = 5, units = "in", res = 600)

scatterplot3d(pc$x[,1:3], pch = pch, color = colors, angle = 20) #log(abs(pc$x[,1:3]))

legend("right", legend = unique(days),
       pch = unique(pch),
       cex = 0.90, 
       horiz = F, inset = -0.0001, xjust = 1)

legend("top", legend = c("T.monococcum","Arabidopsis thaliana", "Zea mays"),
       col =  c("#999999", 
                "#E69F00",
                "#B2182B"),
       cex = 0.90,
       pch = 16,
       inset = -0.0001,
       horiz = T)
dev.off()
library("factoextra")
eig.val <- get_eigenvalue(pc)
var <- get_pca_var(pc)
library("factoextra")
library(ggpubr)
library(gridExtra)

p1 <- fviz_eig(pc, addlabels = TRUE, ylim = c(0, 50))

library("factoextra")
p2 <- fviz_pca_ind(pc, col.ind = "cos2", 
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   repel = TRUE) # Avoid text overlapping (slow if many points)
p2
tiff(filename = "eigen_ind_factoextra.tiff", width = 7, height = 6, units = "in", res = 600)
p2
dev.off()
gridExtra::grid.arrange(p1,p2)

p3 <- fviz_pca_var(pc, col.var = "contrib",
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = T)

