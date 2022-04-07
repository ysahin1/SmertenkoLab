library(magrittr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(DESeq2)

install.packages("DaMiRseq")
source("http://bioconductor.org/biocLite.R")
BiocInstaller::biocLite("vsn")

BiocManager::install(c("HTSFilter"))
source("https://bioconductor.org/biocLite.R")
biocLite("KEGGSOAP")
devtools::install_github(
  c("guangchuangyu/enrichplot",
    "guangchuangyu/DOSE",
    "guangchuangyu/clusterProfiler",
    "guangchuangyu/ChIPseeker"))
devtools::install_github("GuangchuangYu/clusterProfiler")
list.files("../")




###################################################
## load file
###################################################
SRP102142_raw_counts <- read.table("./SRP102142_raw_counts_hisat-sorted.txt", header = TRUE, sep="\t")
head(SRP102142_raw_counts)
row.names(SRP102142_raw_counts) <- SRP102142_raw_counts$Geneid
##########retrieve gene length for GO enrichment
gene_length_SRP102142 <- SRP102142_raw_counts$Length
names(gene_length_SRP102142) <- SRP102142_raw_counts$Geneid
##################################################
SRP102142_raw_counts <- SRP102142_raw_counts[,-c(1:6)] #exclude all information that do not include counts
names(SRP102142_raw_counts) <- c("SRR5358805",
                        "SRR5358806",
                        "SRR5358809",
                        "SRR5358810")
SRP102142_grow_CG <- c("control","control","drought","drought")
SRP102142_grow_CG <- factor(SRP102142_grow_CG, levels = c("control","drought"))
SRP102142_grow_CG <- relevel(SRP102142_grow_CG, ref = "control")
####################################
### generate times series deseq object
####################################
SRP102142_sample_info <- data.frame(condition=SRP102142_grow_CG,
                          row.names=names(SRP102142_raw_counts))

dds_SRP102142 <- DESeqDataSetFromMatrix(countData = SRP102142_raw_counts,
                                   colData = SRP102142_sample_info,
                                   design = ~ condition)

#colData(DESeq.ds)
#assay(DESeq.ds) %>% head

##################################
### filter out count datas and normalization
################################## 
dds_SRP102142 <- dds_SRP102142[rowSums(counts(dds_SRP102142))> 100,]#eleminate counts that are less than zero
length(rowRanges(dds_SRP102142)) # number of genes counted
head(assay(dds_SRP102142))

######normalization
dds_SRP102142 <- estimateSizeFactors(dds_SRP102142)
sizeFactors(dds_SRP102142)
sfn_dds_SRP102142 <- counts(dds_SRP102142, normalized = TRUE)

########transformation
ln_dds_SRP102142 <- log2(sfn_dds_SRP102142+1)


### visulaise norm and log norm data#####
par(mar=c(9,5,3,3))
boxplot(sfn_dds_SRP102142, notch = TRUE,
        main = "Size-Factor Normalized Gene Abundance", 
        ylab = " ",
        las = 2)
axis(1, las=2, outer = FALSE)
par(mar=c(9,5,3,2))
boxplot(ln_dds_SRP102142 , notch = TRUE ,
          main = "Log2-Transformed SizeFactor Normalized Gene Abundance" ,
          ylab = " ",
        las=2)
axis(1, las=2, outer = FALSE)
################# check if your data homoskedastic
plot (ln_dds_SRP102142 [ ,1:2] , cex =.1 , main = " Normalized log2 ( read counts ) " )
################# vst transformation ##########
library(vsn)
vst_dds_SRP102142 <- vst(dds_SRP102142, blind = TRUE)
vstc_dds_SRP102142 <- assay(vst_dds_SRP102142)



######visualize vst norm counts#####
msd_plot <- meanSdPlot (ln_dds_SRP102142,
                           ranks = FALSE, # show the data on the original scale
                           plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing Depth Normalized log2 Transformed Read Counts") +
  ylab("Standard Deviation")
msd_plot <- meanSdPlot(vst_dds_SRP102142,
                        ranks = FALSE, # show the data on the original scale
                        plot = FALSE)
msd_plot$gg +
  ggtitle ( "Estimation of Dispersion Trend") +
  ylab ( "Standard Deviation")

######## hierarchical clustering##########

dis_vstc_dds_SRP102142 <- as.dist(1- cor(vstc_dds_SRP102142, method ="pearson"))

#######visiulaze HC ####
plot(hclust(dis_vstc_dds_SRP102142),
      labels = colnames(dis_vstc_dds_SRP102142),
      xlab = "",
      main = "Hierarchical Clustering of SRP102142")

########### PCA ##########################
pc <- prcomp (t(vstc_dds_SRP102142))

plot(pc$x[ ,1],pc$x[ ,2],
         col = colData(dds_SRP102142)[ ,1],
         main = "Principal Component Analysis")


P <- plotPCA(vst_dds_SRP102142)
############## plot cosmetics ##########################
P <- P + theme_bw() + ggtitle ("Principal Component Analysis of SRP102142")
print (P)

############## DEG analysis ############################
dds_SRP102142 <- DESeq(dds_SRP102142)
dds_SRP102142 <- estimateSizeFactors(dds_SRP102142)
dds_SRP102142 <- estimateDispersions(dds_SRP102142)
dds_SRP102142 <- nbinomWaldTest(dds_SRP102142)

dds_SRP102142_results <- results(dds_SRP102142, 
                                independentFiltering = TRUE, 
                                alpha = 0.05, lfcThreshold = 0)
summary(dds_SRP102142_results)

dn_SRP102142 <- rownames(subset(dds_SRP102142_results, 
                             padj < 0.05 & abs(log2FoldChange) > 0))
length(dn_SRP102142)


# 6 vs 9 degs
#resultsNames(DESeq.ds)

library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(
    A = dn_3vs6,
    B = dn_6vs9,
    C = dn_9vs12,
    D = dn_12vs15,
    E = dn_15vs18
  ),
  filename = "venn_diagramm_overlapDEGs.tiff",
  col = "black",
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  alpha = 0.50,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
          1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.cex = 1.5,
  cat.fontface = "bold",
  margin = 0.05)

library(dplyr)
find.uniques <- list(a=one2two,
                       b=one2three,
                       c=one2four,
                       d=one2five,
                       e=one2six)
as <- find.uniques$a[!(find.uniques$a %in% c(find.uniques$b,
                                   find.uniques$c,
                                   find.uniques$d,
                                   find.uniques$e))]
bs <- find.uniques$b[!(find.uniques$b %in% c(find.uniques$a,
                                             find.uniques$c,
                                             find.uniques$d,
                                             find.uniques$e))]
cs <- find.uniques$c[!(find.uniques$c %in% c(find.uniques$a,
                                             find.uniques$b,
                                             find.uniques$d,
                                             find.uniques$e))]
ds <- find.uniques$d[!(find.uniques$d %in% c(find.uniques$a,
                                             find.uniques$b,
                                             find.uniques$c,
                                             find.uniques$e))]
es <- find.uniques$e[!(find.uniques$e %in% c(find.uniques$a,
                                             find.uniques$b,
                                             find.uniques$c,
                                             find.uniques$d))]


####
library(pheatmap)
topGenes <- head(order(DESeq.ds.results$padj < 0.05),9000)
betas <- coef(DESeq.ds)
betas
mat <- betas[topGenes,-1]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=TRUE, show_rownames = FALSE)

hm.mat_DGEgenes <- log.norm.counts[DEG.names, ]
library(NMF)
aheatmap(hm.mat_DGEgenes,
           Rowv = TRUE, Colv = TRUE,
           distfun = "euclidean", hclustfun = "average" ,
           scale = "row")

par(mar=c(3,4,4,3))
###3 visualize DEGs####
plotMA(DESeq.ds.results, alpha = 0.05 , main = "MA plot of DEGs",
         ylim = c (-8,8))
#write to rownames of counted genes that more then 100 counts
write.table(names(rowRanges(DESeq.ds)), file ="genes_used(more100).txt",
                sep = "\t", quote = FALSE, row.names = FALSE)

######################GO enrichment Analysis (GAGE)####
length(rownames(subset(DESeq.ds.results, padj < 0.05)))
multipleFilters <- readcounts[rownames(subset(DESeq.ds.results, padj < 0.05)),]
head(multipleFilters)
length(multipleFilters$Geneid)
AnnotationOfDEGs.results <- rownames(subset(DESeq.ds.results, padj < 0.05))
head(AnnotationOfDEGs.results)
summary(AnnotationOfDEGs.results)
listDatasets(useMart("ENSEMBL_MART_ENSEMBL")) # list dataset
ensembl <- useMart(biomart = "plants_mart",
                   dataset = "taestivum_eg_gene",
                   host = "plants.ensembl.org")
listFilters(ensembl)
listAttributes(ensembl)

Arabidopsis_homologs <- getBM(attributes = c("ensembl_gene_id",
                                   "athaliana_eg_homolog_ensembl_peptide"),
                    filters = c("ensembl_gene_id"),
                    values = multipleFilters$Geneid,
                    uniqueRows = T,
                    mart = ensembl)# retrive Arabidopsis homologs for 9349 gene_id
summary(Arabidopsis_homologs)
head(Arabidopsis_homologs)
dim(Arabidopsis_homologs)
sum(is.na(Arabidopsis_homologs$athaliana_eg_homolog_ensembl_gene)) # it seems zero
table(Arabidopsis_homologs$athaliana_eg_homolog_ensembl_gene)
Arabidopsis_homologs <- as.data.frame(Arabidopsis_homologs)
Arabidopsis_homologs <- Arabidopsis_homologs[Arabidopsis_homologs$athaliana_eg_homolog_ensembl_gene != '',]
GAGE_DEGs <- rlog.norm.counts[rownames(subset(DESeq.ds.results, padj < 0.05)),] # select DEGs for GAGE object
head(GAGE_DEGs)
anno.DEgenes <- merge(as.data.frame(DESeq.ds.results), features_5, 
                      by.x = "row.names", by.y = "ensembl_gene_id")
anno.DEgenes <- anno.DEgenes[!duplicated(anno.DEgenes$athaliana_eg_homolog_ensembl_gene),]
anno.DEgenes <- merge(rlog.norm.counts, anno.DEgenes, by.x="row.names", by.y ="Row.names")
row.names(anno.DEgenes) <- anno.DEgenes$athaliana_eg_homolog_ensembl_gene
anno.DEgenes <- anno.DEgenes[,-1]
anno.DEgenes <- anno.DEgenes[,-c(31,32,33,34,35,36,37)]
#####################AGRİGO#################################
length(rownames(subset(DESeq.ds.results , padj < 0.05))) #number of DEGs
agri.go <- rownames(subset(DESeq.ds.results , padj < 0.05))
library("biomaRt") #Annotation
listDatasets(useMart("ENSEMBL_MART_ENSEMBL")) # list dataset
ensembl <- useMart(biomart = "plants_mart",
                   dataset = "taestivum_eg_gene",
                   host = "plants.ensembl.org")
features.agriGo <- getBM(attributes = c("ensembl_gene_id",
                                 "go_id"),
                  filters = c("ensembl_gene_id"),
                  values = agri.go, 
                  mart = ensembl)

features.agriGo <- features.agriGo[features.agriGo$go_id != "",]

write.table(features.agriGo, file = "agrigo.txt",sep = "\t", quote = FALSE, row.names = FALSE)
######################################################
####ANNOTATION#####

library("biomaRt") #Annotation
listDatasets(useMart("ENSEMBL_MART_ENSEMBL")) # list dataset
ensembl <- useMart(biomart = "plants_mart",
                   dataset = "taestivum_eg_gene",
                   host = "plants.ensembl.org")
listFilters(ensembl) # annotation attributes
DGEgenes <- rownames(DESeq.ds.results) # for biomart

features <- getBM(attributes = c("ensembl_gene_id",
                                 "go_id"),
                  filters = c("ensembl_gene_id"),
                  values = DGEgenes, 
                  mart = ensembl)

##################### TOPGO enrichment analysis################################
library(topGO) 
features <- as.data.frame(features)
features <- features[features$go_id !="",]

geneID2GO <- by(features$go_id,
                features$ensembl_gene_id,
                function(x) as.character(x))

all.genes <- sort(unique(as.character(features$ensembl_gene_id)))
int.genes <-  rownames(subset(DESeq.ds.results, padj < 0.05))# sig. genes
int.genes <- factor(as.integer(all.genes %in% int.genes))
names(int.genes) = all.genes

go.obj <- new("topGOdata", ontology='BP'
              , allGenes = int.genes
              , nodeSize = 10
              , annot = annFUN.gene2GO
              , gene2GO = geneID2GO)

resultFisher <- runTest(go.obj, algorithm = "classic", statistic = "fisher")

# eleminate less than 0.05
sig.go <- score(resultFisher)[score(resultFisher)< 0.05] 
length(sig.go)
?par 
par(cex = 0.05)
GO.graph <-  showSigOfNodes(go.obj, sig.go, firstSigNodes = 10, useInfo = 'all', sigForAll = TRUE, swPlot = TRUE)
GO.graph$dag


go.result <- GenTable(go.obj, classicFisher = resultFisher, ranksOf = "classicFisher")
go.result
goID <- go.result[2,"GO.ID"]

print(showGroupDensity(go.obj, goID, ranks = TRUE, rm.one = TRUE))

printGenes(go.obj, whichTerms = goID, file="hebelehubele.txt") # çalışmıyor


############ Pathway /cluster enrichment analysis####################
listAttributes(ensembl)
DEgenes_FC <- rownames(DESeq.ds.results)
features_5 <- getBM(attributes = c("ensembl_gene_id",
                                   "athaliana_eg_homolog_ensembl_gene"),
                    filters = c("ensembl_gene_id"),
                    values = DEgenes_FC, 
                    mart = ensembl)
features_5 <- as.data.frame(features_5)
features_5 <- features_5[features_5$athaliana_eg_homolog_ensembl_gene != '',]

anno <-  select(org.At.tair.db, keys = as.vector(features_5$athaliana_eg_homolog_ensembl_gene), 
                keytype = "TAIR", columns = c("ENTREZID"))
keytypes(org.At.tair.db)
anno <- as.data.frame(anno)
anno <- anno[anno$ENTREZID != '',]
# merge data to our matrix
out.gene.FC <- merge(as.data.frame(DESeq.ds.results), features_5, 
                     by.x = "row.names", by.y = "ensembl_gene_id")
out.gene.FC <- merge(out.gene.FC, anno, by.x="athaliana_eg_homolog_ensembl_gene", by.y ="TAIR")


out.gene.FC <- out.gene.FC[out.gene.FC$ENTREZID != '',]
out.gene.FC <- out.gene.FC[!duplicated(out.gene.FC$ENTREZID),]# drop null data
##### Worked Pathwayview ##
out.gene.FC <- na.omit(out.gene.FC)
out.gene.FC
out.gene.FC.clust <- out.gene.FC[abs(out.gene.FC$log2FoldChange) > 2 & out.gene.FC$padj < 0.05,]
out.gene.FC.clust_2 <- out.gene.FC.clust$padj
names(out.gene.FC.clust_2) <- out.gene.FC.clust$ENTREZID 
out.gene.FC.clust_2 <- na.exclude(out.gene.FC.clust_2)
sorted.genes <- sort(out.gene.FC.clust_2, decreasing = T)
library(pathview)
pv.out <- pathview(sorted.genes, pathway.id = "04146",
                   species = "ath", out.suffix = "amk-peroxisome", kegg.native = T,
                   same.layer = F)
warning(pv.out)

head(pv.out$plot.data.gene)
anno[as.vector(grep("ALDH2C",anno$SYMBOL )),]

BiocManager::install("enrichplot")
library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(pathview)
sorted.genes <- sort(out.gene.FC.clust_2, decreasing = TRUE)
de <- names(sorted.genes)[abs(sorted.genes) > 6]
ego <- enrichGO(de, OrgDb = "org.At.tair.db", ont="BP", readable=TRUE)
ego
goplot(ego, showCategory = 5)
?goplot
barplot(ego, showCategory=10)
ego2 <- simplify(ego)
cnetplot(ego2, foldChange=sorted.genes, circular = TRUE, colorEdge = TRUE, showCategory = 5)
?cnetplot
go <- enrichGO(de, OrgDb = "org.At.tair.db", ont="ALL")
dotplot(de, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
heatplot(ego2, foldChange=sorted.genes, showCategory = 10)
?heatplot
emapplot(ego2)

library("biomaRt") #Annotation
listDatasets(useMart("ENSEMBL_MART_ENSEMBL")) # list dataset
ensembl <- useMart(biomart = "plants_mart",
                   dataset = "taestivum_eg_gene",
                   host = "plants.ensembl.org")
listFilters(ensembl) # annotation attributes

head(read.counts_2)
anno_reads <- getBM(attributes = c("ensembl_gene_id",
                                   "athaliana_eg_homolog_ensembl_gene"),
                    filters = c("ensembl_gene_id"),
                    values = rownames(read.counts_2), 
                    mart = ensembl)
#############DATA MUNGING#####
table(is.na(features$uniprotswissprot))
features_1 <- features[features$uniprotswissprot != '',] # remove blank entries
table(is.na(features_1))

out.gene.FC <- merge(as.data.frame(DESeq.ds.results), features_1, 
                     by.x = "row.names", by.y = "ensembl_gene_id")

out.gene.FC <- as.data.frame(out.gene.FC)
dim(out.gene.FC)

out.gene.FC <- out.gene.FC[out.gene.FC$uniprotswissprot != '',]
dim(out.gene.FC)

hebele <- cbind(out.gene.FC$log2FoldChange,out.gene.FC$uniprotswissprot)
row.names(hebele) <- hebele[,2]
colnames(hebele) <- c("log2FoldChange", "uniprotswissprot")
hebele <- hebele[,1]
length(hebele)
write.table(rownames(hebele), file = "uniportnames.txt",sep = "/t", quote = FALSE, row.names = FALSE)
convertId(hebele, dataset = "taestivum_eg_gene", filters = "uniprotswissprot", attributes =C(filters, "entrezgene_id"), keepNoId = T, keepMultipleId = F, verbose = F)

pv.out <- pathview(hebele, pathway.id = "00940",
                   species = "osa", out.suffix = "at_deneme", kegg.native = T, same.layer = F)
head(pv.out$plot.data.gene)
list.files(pattern="osa00190", full.names=T)
str(pv.out)
head(row.names(subset(hebele, log2FoldChange >= 1)))

download_KEGGfile(pathway_id="01100", species= 'ath')

genes<-rownames(hebele)
genes <- na.omit(anno$ENTREZID)
?find_enriched_pathway

pho_KEGGresult<-find_enriched_pathway(hebele, species='ath',  download_latest = TRUE)
pho_KEGGresult
pho_KEGGresult[[1]][,c(1,5)]

columns(KEGG.db)





#####Single gene visilualize#########
write.table(features$ensembl_gene_id, file ="ensembl_gene_id.txt",
            sep = ",", quote = FALSE, row.names = FALSE)

des_fea <- out.gene.FC.clust$GENENAME
grep("3-keto", des_fea)
out.gene.FC.clust[782,]
AGO1 <- "TraesCS7A02G040700"
AGO4 <- "TraesCS3A02G188400"
Actin <- "TraesCS4A02G334700"
EFa <- "TraesCS2A02G083300"
EFa1 <- "TraesCS4A02G107700"
GG <- "TraesCS4A02G175900"
CDCP <- "TraesCS3A02G369000"
AcS <- "TraesCS2A02G290900"
keto <- "TraesCS4A02G007400"
mito <- "TraesCS7A02G335300"
pect <- "TraesCS2A02G109600"
eps <- "TraesCS3A02G118400"

plotCounts ( dds = DESeq.ds ,
             gene = AGO1 ,
             normalized = TRUE , transform = FALSE ,
             main = expression(atop("Expression of AGO1")))
plotCounts ( dds = DESeq.ds ,
             gene = AGO4 ,
             normalized = TRUE , transform = FALSE ,
             main = expression(atop("Expression of AGO4")))
plotCounts ( dds = DESeq.ds ,
             gene = CDCP ,
             normalized = TRUE , transform = FALSE ,
             main = expression(atop("Expression of CDCP")))
plotCounts ( dds = DESeq.ds ,
             gene = Actin ,
             normalized = TRUE , transform = FALSE ,
             main = expression(atop("Expression of Actin")))
plotCounts ( dds = DESeq.ds ,
             gene = EFa1 ,
             normalized = TRUE , transform = FALSE ,
             main = expression(atop("Expression of EFa1")))
plotCounts ( dds = DESeq.ds ,
             gene = GG ,
             normalized = TRUE , transform = FALSE ,
             main = expression(atop("Expression of Glucan")))
plotCounts ( dds = DESeq.ds ,
             gene = AcS ,
             normalized = TRUE , transform = FALSE ,
             main = expression(atop("Expression of AcS")))
plotCounts ( dds = DESeq.ds ,
             gene = keto ,
             normalized = TRUE , transform = FALSE ,
             main = expression(atop("3-ketoacyl-CoA synthase 6")))
plotCounts ( dds = DESeq.ds ,
             gene = mito ,
             normalized = TRUE , transform = FALSE ,
             main = expression(atop("mitogen")))
plotCounts ( dds = DESeq.ds ,
             gene = eps ,
             normalized = TRUE , transform = FALSE ,
             main = expression(atop("eps")))
DGE.results.sorted_logFC <- DESeq.ds.results[order(DESeq.ds.results$log2FoldChange),]
DGEgenes_logFC <- rownames(subset(DGE.results.sorted_logFC,padj < 0.05))
head(features[match(DGEgenes_logFC, features[16123:16162,2] ),])
      