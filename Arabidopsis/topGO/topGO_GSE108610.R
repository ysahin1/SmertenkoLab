library(topGO)
GO_slim_table <- read.delim2("/media/yunus/TOSHIBA1TB/c3-c4/peroxisome-C3vsC4/DEG/arabidopsis/ATH_GO_GOSLIM.txt", header=FALSE, comment.char="#")
head(GO_slim_table)
GOs_GSE108610<- data.frame(gene_id=GO_slim_table$V1,
                          go=GO_slim_table$V6)
geneID2GO <- by(GOs_GSE108610$go,
                GOs_GSE108610$gene_id,
                function(x) as.character(x))
# create named factor vector
all.genes <- sort(unique(as.character(GOs_GSE108610$gene_id)))
int.genes <-  rownames(subset(dds_GSE108610_results, padj < 0.05))# sig. genes
int.genes <- factor(as.integer(all.genes %in% int.genes))
names(int.genes) = all.genes
#####3 create topGo object

############ BP ######################
go_BP_GSE108610<- new("topGOdata"
                     , ontology='BP'
                     , allGenes = int.genes
                     , nodeSize = 1
                     , annot = annFUN.gene2GO
                     , gene2GO = geneID2GO)

### make your analysis
resultFisher_BP_GSE108610<- runTest(go_BP_GSE108610, 
                                   algorithm = "classic", 
                                   statistic = "fisher")
number_bp_top_nodes <- length(score(resultFisher_BP_GSE108610))
result_BP_GSE108610<- GenTable(go_BP_GSE108610, 
                              classicFisher = resultFisher_BP_GSE108610, 
                              orderBy = "classicFisher", 
                              ranksOf = "classicFisher", 
                              topNodes = number_bp_top_nodes)
length(score(resultFisher_BP_GSE108610))
result_BP_GSE108610$adj.pvalue <- p.adjust(score(resultFisher_BP_GSE108610),
                                          method="BH")
result_BP_GSE108610$ontology <- rep("BP", 
                                   dim(result_BP_GSE108610)[1])
result_BP_GSE108610$ratio <- result_BP_GSE108610$Significant/result_BP_GSE108610$Annotated

# select those terms with a p-value < 0.05
result_BP_GSE108610$genes <- sapply(result_BP_GSE108610$GO.ID, function(x)
{
  genes<-genesInTerm(go_BP_GSE108610, x)
  genes[[1]][genes[[1]] %in% dn_GSE108610] # myGenes is the queried gene list
})
head(result_BP_GSE108610)
enriched_BP_GSE108610<- result_BP_GSE108610[result_BP_GSE108610$adj.pvalue < .05,]
#Number of Enriched BP_GO terms
dim(enriched_BP_GSE108610)[1]
write.table(enriched_BP_GSE108610[,-c(10)], file = "/home/yunus/Masaüstü/Mega/peroxisome-C3vsC4/DEG/arabidopsis/for_article/topgo/tables/rice/GSE108610_topgo_BP_res.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
# Bu çalışmada GO enrichment açısından önemli olan gen sayısı
table(dn_GSE108610%in% unique(unlist(result_BP_GSE108610$genes)))[2]
################# MF #####################
go_MF_GSE108610<- new("topGOdata"
                     , ontology='MF'
                     , allGenes = int.genes
                     , nodeSize = 1
                     , annot = annFUN.gene2GO
                     , gene2GO = geneID2GO)

### make your analysis
resultFisher_MF_GSE108610<- runTest(go_MF_GSE108610, 
                                   algorithm = "classic", 
                                   statistic = "fisher")
number_MF_top_nodes <- length(score(resultFisher_MF_GSE108610))
result_MF_GSE108610<- GenTable(go_MF_GSE108610, 
                              classicFisher = resultFisher_MF_GSE108610, 
                              orderBy = "classicFisher", 
                              ranksOf = "classicFisher", 
                              topNodes = number_MF_top_nodes)
length(score(resultFisher_MF_GSE108610))
result_MF_GSE108610$adj.pvalue <- p.adjust(score(resultFisher_MF_GSE108610),
                                          method="BH")
result_MF_GSE108610$ontology <- rep("MF", 
                                   dim(result_MF_GSE108610)[1])
result_MF_GSE108610$ratio <- result_MF_GSE108610$Significant/result_MF_GSE108610$Annotated

# select those terms with a p-value < 0.05
result_MF_GSE108610$genes <- sapply(result_MF_GSE108610$GO.ID, function(x)
{
  genes<-genesInTerm(go_MF_GSE108610, x)
  genes[[1]][genes[[1]] %in% dn_GSE108610] # myGenes is the queried gene list
})
head(result_MF_GSE108610)
enriched_MF_GSE108610<- result_MF_GSE108610[result_MF_GSE108610$adj.pvalue < .05,]
#Number of Enriched BP_GO terms
dim(enriched_MF_GSE108610)[1]
# Bu çalışmada GO enrichment açısından önemli olan gen sayısı
table(dn_GSE108610%in% unique(unlist(result_MF_GSE108610$genes)))[2]
length(enriched_MF_GSE108610$genes[["X"]]) # ne istersen bal GO term yaz.
enriched_MF_GSE108610[,1:8]
#######################CC############################################

go_CC_GSE108610<- new("topGOdata"
                     , ontology='CC'
                     , allGenes = int.genes
                     , nodeSize = 1
                     , annot = annFUN.gene2GO
                     , gene2GO = geneID2GO)

### make your analysis
resultFisher_CC_GSE108610<- runTest(go_CC_GSE108610, 
                                   algorithm = "classic", 
                                   statistic = "fisher")
number_CC_top_nodes <- length(score(resultFisher_CC_GSE108610))
result_CC_GSE108610<- GenTable(go_CC_GSE108610, 
                              classicFisher = resultFisher_CC_GSE108610, 
                              orderBy = "classicFisher", 
                              ranksOf = "classicFisher", 
                              topNodes = number_CC_top_nodes)
length(score(resultFisher_CC_GSE108610))
result_CC_GSE108610$adj.pvalue <- p.adjust(score(resultFisher_CC_GSE108610),
                                          method="BH")
result_CC_GSE108610$ontology <- rep("CC", 
                                   dim(result_CC_GSE108610)[1])
result_CC_GSE108610$ratio <- result_CC_GSE108610$Significant/result_CC_GSE108610$Annotated

# select those terms with a p-value < 0.05
result_CC_GSE108610$genes <- sapply(result_CC_GSE108610$GO.ID, function(x)
{
  genes<-genesInTerm(go_CC_GSE108610, x)
  genes[[1]][genes[[1]] %in% dn_GSE108610] # myGenes is the queried gene list
})
head(result_CC_GSE108610)
enriched_CC_GSE108610<- result_CC_GSE108610[result_CC_GSE108610$adj.pvalue < .05,]
#Number of Enriched BP_GO terms
dim(enriched_CC_GSE108610)[1]
# Bu çalışmada GO enrichment açısından önemli olan gen sayısı
table(dn_GSE108610%in% unique(unlist(result_CC_GSE108610$genes)))[2]
length(enriched_CC_GSE108610$genes[["X"]]) # ne istersen bal GO term yaz.
enriched_CC_GSE108610[,1:8]







enriched_TOPGOs_GSE108610<- rbind(enriched_BP_GSE108610, 
                                 enriched_MF_GSE108610, 
                                 enriched_CC_GSE108610)
dim(enriched_TOPGOs_GSE108610)
####################visualize results###################
########################################################
names(enriched_TOPGOs_GSE108610) <- c("category",
                                     "term",
                                     "numInCat",
                                     "Gene_Number",
                                     "Expected",
                                     "over_represented_pvalue",
                                     "Q_value",
                                     "ontology",
                                     "Rich_Factor",
                                     "Genes")
dottplot <- function(df, showCategory){
  df <- df[with(df, order(Rich_Factor, Q_value, decreasing = c(TRUE, FALSE))),]
  df <- head(df, n=showCategory)
  d_plot <- ggplot(df, aes_string(x="term", 
                                  y="Rich_Factor", 
                                  colour="Q_value",
                                  size="Gene_Number")) + 
    geom_point() +
    facet_grid(~ontology)+
    scale_color_gradient(low="#FF0000",
                         high="#000000") +
    coord_flip() +
    theme_bw(base_size=9)
  return(d_plot)
} 
library(ggplot2)
dottplot(enriched_TOPGOs_GSE108610, 
         showCategory = 107)



library(GOplot)
DEGs.table <- dds_GSE108610_results[dn_GSE108610,]
genes=DEGs.table
terms=enriched_TOPGOs_GSE108610
names(terms) <- c("Category",
                  "term",
                  "numInCat",
                  "numDEInCat",
                  "under_represented_pvalue",
                  "over_represented_pvalue",
                  "adj_pval",
                  "ontology",
                  "Rich_Factor",
                  "Genes"
                  
)
names(genes) <- c("baseMean",
                  "logFC",
                  "lfcSE",
                  "stat",
                  "P.Value",
                  "adj.P.Val")
genes$ID <- rownames(genes)
logFC <- sapply(as.vector(unlist(terms$Genes)), function(x) genes$logFC[match(x, genes$ID)])
count <- sapply(1:length(names(terms$Genes)), function(x) length(terms$Genes[[x]]))
#logFC[is.na(logFC)] <- 0 #Fold change değeri olmayanlara sıfır değeri verdim 

s <- 1; zsc <- c()
for (c in 1:length(count)){
  value <- 0
  e <- s + count[c] - 1
  value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 1, -1))
  zsc <- c(zsc, sum(value)/sqrt(count[c]))
  s <- e + 1
}

df_for_go_circle <- data.frame(category = rep(as.character(terms$Category), count), 
                               term = rep(as.character(terms$term), count),
                               genes = as.vector(unlist(terms$Genes)), 
                               logFC = as.vector(logFC), 
                               adj_pval = rep(terms$adj_pval, count),
                               zscore = rep(zsc, count), 
                               stringsAsFactors = FALSE)



enriched_BP_GSE108610[enriched_BP_GSE108610$Term == "regulation of secondary metabolite biosy...",-10]
unique(df_for_go_circle$term[order(df_for_go_circle$term)])
round(min(df_for_go_circle$logFC[df_for_go_circle$term == "regulation of secondary metabolite biosy..."]), digits = 2)
round(max(df_for_go_circle$logFC[df_for_go_circle$term == "regulation of secondary metabolite biosy..."]), digits = 2)


