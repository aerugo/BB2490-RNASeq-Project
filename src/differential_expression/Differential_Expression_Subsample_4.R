### Unix codes to merge files together
## join -j 1  SN_10_LPS_CGATGT_L001_count.txt SN_10_LPS_CGATGT_L002_count.txt| join -j 1 SN_10_LPS_CGATGT_L003_count.txt -  | join -j 1 SN_10_LPS_CGATGT_L004_count.txt - > SN10_LPS_Subsample4.txt
## join -j 1  SN10_UNST_ATCACG_L001_count.txt SN10_UNST_ATCACG_L002_count.txt| join -j 1 SN10_UNST_ATCACG_L003_count.txt -  | join -j 1 SN10_UNST_ATCACG_L004_count.txt - > SN10_UNST_Subsample4.txt
## join -j 1  SN11_UNST_TTAGGC_L002_count.txt SN11_UNST_TTAGGC_L003_count.txt| join -j 1 SN11_UNST_TTAGGC_L004_count.txt -  | join -j 1 SN11_UNST_TTAGGC_L005_count.txt - > SN11_UNST_Subsample4.txt
## join -j 1  SN_11_LPS_TGACCA_L002_count.txt SN_11_LPS_TGACCA_L003_count.txt| join -j 1 SN_11_LPS_TGACCA_L004_count.txt -  | join -j 1 SN_11_LPS_TGACCA_L005_count.txt - > SN11_LPS_Subsample4.txt
## join -j 1  SN12_UNST_ACAGTG_L003_count.txt SN12_UNST_ACAGTG_L004_count.txt| join -j 1 SN12_UNST_ACAGTG_L005_count.txt -  | join -j 1 SN12_UNST_ACAGTG_L006_count.txt - > SN12_UNST_Subsample4.txt
## join -j 1  SN_12_LPS_GCCAAT_L003_count.txt SN_12_LPS_GCCAAT_L004_count.txt| join -j 1 SN_12_LPS_GCCAAT_L005_count.txt -  | join -j 1 SN_12_LPS_GCCAAT_L006_count.txt - > SN12_LPS_Subsample4.txt

#######################################################################
## Differential expresssion concencerning 100% of data from flowcells
## DEQSeq2 : Differential gene expression analysis in Inflammatoy and Non Inflammatory agents in WBC ##

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")
## Differential expression analysis for Sample S10 ##ß

setwd("~/BB2490-RNASeq-Project/data/Htseq")
SN10_LPS_Subsample4 <- read.table("~/BB2490-RNASeq-Project/data/Htseq/SN10_LPS_Subsample4.txt", quote="\"", comment.char="")
SN10_UNST_Subsample4 <- read.table("~/BB2490-RNASeq-Project/data/Htseq/SN10_UNST_Subsample4.txt", quote="\"", comment.char="")

colnames (SN10_UNST_Subsample4) <- c("Gene_features", "UNST_Lane_1", "UNST_Lane_2", "UNST_Lane_3", "UNST_Lane_4")
colnames (SN10_LPS_Subsample4) <- c("Gene_features", "LPS_Lane_1", "LPS_Lane_2", "LPS_Lane_3", "LPS_Lane_4") 
data_SN10 <- merge (SN10_LPS_Subsample4, SN10_UNST_Subsample4, by.y ="Gene_features" )
## Considering each lanes as technical replicates and hence given names based on lane numbers

N = dim(data_SN10)[1]
rownames(data_SN10) = data_SN10[,1]
data_SN10 = data_SN10[,-1]
data_SN10= data_SN10[c(6:N),]
## removing last 5 rows which in our case turn out be in top 5 rows
data_SN10 <- data_SN10[ rowSums(data_SN10) > 1, ] 
## Filtering to reduce number of genes that have 0 count values
## 25347 ENSEMBL genes

colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS", "LPS","UNST", "UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data_SN10, colData, formula(~ condition))
dds <- DESeq(dds)
# plotMA(dds, main="Differential Gene Expression in Sample S10 at Subsample 100% data")


pdf("../../results/MA_Plot_SN10_SB4.pdf", height=8, width=12)
plotMA(dds, main="Differential Gene Expression in Sample S10 at Subsample 100% data")
dev.off()
## The MA plot kinda of supports the argument of large number of genes differential expressed between two condition

res <- results(dds)
res_clean <- res[(!is.na(res$padj)) &
                        (res$padj != 0.000000e+00), ]
## I did this filtering to remove genes with 0 padjusted values
## May be it would be interesting to see why there is padjusted to 0

resOrdered <- res_clean[order(res_clean$padj),]
head(resOrdered)
sig <- resOrdered[resOrdered$padj<0.05 &
                    abs(resOrdered$log2FoldChange)>=1,]
## This gave us 1804 genes that are differential expressed with the above 
## criteria.
## Defning the criteria for the genes which are significant as ones
## with the Padjusted values lesser than 5% FDR and havin log2Fold change
## greater than 1. It would be interesting to see what happens with different
## cutoff. The above results ga

summary(res_clean)
#4323 genes upregulated and 4323 genes downregulated
sum(res_clean$padj < 0.1, na.rm=TRUE)
## 8646 genes are differentially expressed at 10% FDR

sum(res_clean$padj < 0.05, na.rm=TRUE)
## 7747 genes are differentially expressed at 1% FDR

S10_genes_sc <- as.character(sig@rownames)


############################################################################################
## Differential expression analysis for Sample S11 ##

setwd("~/BB2490-RNASeq-Project/data/Htseq")
SN11_LPS_Subsample4 <- read.table("~/BB2490-RNASeq-Project/data/Htseq/SN11_LPS_Subsample4.txt", quote="\"", comment.char="")
SN11_UNST_Subsample4 <- read.table("~/BB2490-RNASeq-Project/data/Htseq/SN11_UNST_Subsample4.txt", quote="\"", comment.char="")

colnames (SN11_LPS_Subsample4) <- c("Gene_features", "LPS_Lane_1", "LPS_Lane_2", "LPS_Lane_3", "LPS_Lane_4")
colnames (SN11_UNST_Subsample4) <- c("Gene_features", "UNST_Lane_1", "UNST_Lane_2", "UNST_Lane_3", "UNST_Lane_4" ) 
data_SN11 <- merge (SN11_LPS_Subsample4, SN11_UNST_Subsample4, by.y ="Gene_features"  )

N = dim(data_SN11)[1]
rownames(data_SN11) = data_SN11[,1]
data_SN11 = data_SN11[,-1]
data_SN11 = data_SN11[c(6:N),]
data_SN11 <- data_SN11[ rowSums(data_SN11) > 1, ]

colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS", "LPS","UNST", "UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data_SN11, colData, formula(~ condition))
dds <- DESeq(dds)


pdf("../../results/MA_Plot_SN11_SB4.pdf", height=8, width=12)
plotMA(dds, main="Differential Gene Expression in Sample S11 at Subsample 100% data")
dev.off()

## See the above argument to follow the steps done in the following steps
res_SN11 <- results(dds)
res_clean_SN11 <- res_SN11[(!is.na(res_SN11$padj)) &
                   (res_SN11$padj != 0.000000e+00), ]
resOrdered_SN11 <- res_clean_SN11[order(res_clean_SN11$padj),]
head(resOrdered_SN11)
sig_S11 <- resOrdered_SN11[resOrdered_SN11$padj<0.05 &
                    abs(resOrdered_SN11$log2FoldChange)>=1,]

## This gave us 1169 genes that are differentially expressed using above criteria values



summary(res_clean_SN11)
## 4068 genes upregulated and 4367 genes downregulated
sum(res_clean_SN11$padj < 0.1, na.rm=TRUE)
## 8345 genes differentially expressed
sum(res_clean_SN11$padj < 0.05, na.rm=TRUE)
## 7557 genes differentially expressed

S11_genes_sc <- as.character(sig_S11@rownames)

###########################################################
## Differential expression analysis for Sample S12 ##

setwd("~/BB2490-RNASeq-Project/data/Htseq")
SN12_LPS_Subsample4 <- read.table("~/BB2490-RNASeq-Project/data/Htseq/SN12_LPS_Subsample4.txt", quote="\"", comment.char="")
SN12_UNST_Subsample4 <- read.table("~/BB2490-RNASeq-Project/data/Htseq/SN12_UNST_Subsample4.txt", quote="\"", comment.char="")

colnames (SN12_LPS_Subsample4) <- c("Gene_features", "LPS_Lane_1", "LPS_Lane_2", "LPS_Lane_3", "LPS_Lane_4")
colnames (SN12_UNST_Subsample4) <- c("Gene_features", "UNST_Lane_1", "UNST_Lane_2", "UNST_Lane_3", "UNST_Lane_4" ) 
data_S12 <- merge (SN12_LPS_Subsample4, SN12_UNST_Subsample4, by.y ="Gene_features"  )

N = dim(data_S12)[1]
rownames(data_S12) = data_S12[,1]
data_S12 = data_S12[,-1]
data_S12 = data_S12[c(6:N),]
data_S12 <- data_S12[ rowSums(data_S12) > 1, ]

colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS", "LPS","UNST", "UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data_S12, colData, formula(~ condition))
dds <- DESeq(dds)

pdf("../../results/MA_Plot_SN12_SB4.pdf", height=8, width=12)
plotMA(dds, main="Differential Gene Expression in Sample S12 at Subsample 100% data")
dev.off()

res_S12 <- results(dds)
res_clean_SN12 <- res_S12[(!is.na(res_S12$padj)) &(res_S12$padj != 0.000000e+00), ]
resOrdered_SN12 <- res_clean_SN12[order(res_clean_SN12$padj),]
head(resOrdered_SN12)
sig_S12 <- resOrdered_SN12[resOrdered_SN12$padj<0.05 & abs(resOrdered_SN12$log2FoldChange)>=1,]
## 897 genes are ther that is differentially expressed


summary(res_S12)
## 3643 genes are upregulated and 4099 genes are downregulated
sum(res_S12$padj < 0.1, na.rm=TRUE)
## 7742 genes DExpressed

sum(res_S12$padj < 0.05, na.rm=TRUE)
## 6895 DExpressed

S12_genes_sc <- as.character(sig_S12@rownames)



### Visualizing the common genes number identified in each of samples at different depth

## Finding the common genes that are differntially expressed and havníng fold change>=1


library(gplots)
library(VennDiagram)
Common_genes<- Reduce(intersect,  list(S10_genes_sc, 
                        S11_genes_sc,
                        S12_genes_sc))
## 516 genes 

venn.plot<-venn.diagram(x = list(SN10=S10_genes_sc, 
                                 SN11=S11_genes_sc,
                                 SN12=S12_genes_sc),
                        filename = "../../results/Venn_plot.tiff",
                        col = "black",
                        lty = "blank",
                        lwd = 4,
                        fill = c("cornflowerblue", "green", "yellow"),
                        alpha = 0.50,
                        cex = 1.5,
                        fontfamily = "serif",
                        fontface = "bold",
                        #cat.pos = c(-0,0,-0),
                        cat.col = c("darkblue", "darkgreen", "orange"),
                        cat.cex = 1.0)

