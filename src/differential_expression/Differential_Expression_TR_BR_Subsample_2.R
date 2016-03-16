## Differential Expression analysis for Subsample 2

###### Awk codes for combining files together  for technical replicates (Possible for cooler solutions)
## join -j 1  SN_10_LPS_CGATGT_L002_count.txt SN_10_LPS_CGATGT_L001_count.txt| join -j 1 SN_10_LPS_CGATGT_L003_count.txt -  | join -j 1 SN_10_LPS_CGATGT_L004_count.txt - > SN10_LPS_SS2.txt
## join -j 1  SN10_UNST_ATCACG_L002_count.txt SN10_UNST_ATCACG_L001_count.txt| join -j 1 SN10_UNST_ATCACG_L003_count.txt -  | join -j 1 SN10_UNST_ATCACG_L004_count.txt - > SN10_UNST_SS2.txt
## join -j 1  SN11_UNST_TTAGGC_L003_count.txt SN11_UNST_TTAGGC_L002_count.txt| join -j 1 SN11_UNST_TTAGGC_L004_count.txt -  | join -j 1 SN11_UNST_TTAGGC_L005_count.txt - > SN11_UNST_SS2.txt
## join -j 1  SN_11_LPS_TGACCA_L003_count.txt SN_11_LPS_TGACCA_L002_count.txt| join -j 1 SN_11_LPS_TGACCA_L004_count.txt -  | join -j 1 SN_11_LPS_TGACCA_L005_count.txt - > SN11_LPS_SS2.txt
## join -j 1  SN_12_LPS_GCCAAT_L004_count.txt SN_12_LPS_GCCAAT_L003_count.txt| join -j 1 SN_12_LPS_GCCAAT_L005_count.txt -  | join -j 1 SN_12_LPS_GCCAAT_L006_count.txt - > SN12_LPS_SS2.txt
##join -j 1  SN12_UNST_ACAGTG_L004_count.txt SN12_UNST_ACAGTG_L003_count.txt| join -j 1 SN12_UNST_ACAGTG_L005_count.txt -  | join -j 1 SN12_UNST_ACAGTG_L006_count.txt - > SN12_UNST_SS2.txt

## For Biological Replicates
# awk '{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print $1 " " sum}' SN10_LPS_SS2.txt > LPS_temp1
# awk '{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print $1 " " sum}' SN11_LPS_SS2.txt > LPS_temp2
# awk '{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print $1 " " sum}' SN12_LPS_SS2.txt > LPS_temp3
# 
# awk '{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print $1 " " sum}' SN10_UNST_SS2.txt > UNST_temp1
# awk '{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print $1 " " sum}' SN11_UNST_SS2.txt > UNST_temp2
# awk '{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print $1 " " sum}' SN12_UNST_SS2.txt > UNST_temp3

# join -j 1  LPS_temp2 LPS_temp3| join -j 1 LPS_temp1 -  > LPS_SS2.txt
# join -j 1  UNST_temp2 UNST_temp3| join -j 1 UNST_temp1 -  > UNST_SS2.txt

####################################################################################################################
## Differential expression for technical replicates considering 

library("DESeq2")
####################################################################################################################
## Differential expression analysis for Sample S10 ##

setwd("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_2")
SN10_LPS_SS2 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_2/SN10_LPS_SS2.txt", quote="\"", comment.char="")
SN10_UNST_SS2 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_2/SN10_UNST_SS2.txt", quote="\"", comment.char="")

colnames (SN10_LPS_SS2) <- c("Gene_features", "LPS_Lane_4", "LPS_Lane_3", "LPS_Lane_2", "LPS_Lane_1")
colnames (SN10_UNST_SS2) <- c("Gene_features", "UNST_Lane_4", "UNST_Lane_3", "UNST_Lane_2", "UNST_Lane_1" ) 
data_SN10_SS2 <- merge (SN10_LPS_SS2, SN10_UNST_SS2, by.y ="Gene_features")

N = dim(data_SN10_SS2)[1]
rownames(data_SN10_SS2) = data_SN10_SS2[,1]
data_SN10_SS2 = data_SN10_SS2[,-1]
data_SN10_SS2 = data_SN10_SS2[c(6:N),]
data_SN10_SS2 <- data_SN10_SS2[ rowSums(data_SN10_SS2) > 1,]
## 25,347 genes in consideration

colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS", "LPS","UNST", "UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data_SN10_SS2, colData, formula(~ condition))
dds <- DESeq(dds)


pdf("~/BB2490-RNASeq-Project/results/MA_plots/MA_technical_replicates_SN10_SS2.pdf", height=8, width=12)
plotMA(dds, main="Differential Gene Expression in S10 techincal replicates at Subsample 2")
dev.off()

## See the above argument to follow the steps done in the following steps
res_SN10_SS2 <- results(dds)
res_clean_SN10_SS2 <- res_SN10_SS2[(!is.na(res_SN10_SS2$padj)) &
                             (res_SN10_SS2$padj != 0.000000e+00), ]
resOrdered_SN10_SS2 <- res_clean_SN10_SS2[order(res_clean_SN10_SS2$padj),]
head(resOrdered_SN10_SS2)
sig_S10_SS2 <- resOrdered_SN10_SS2[resOrdered_SN10_SS2$padj<0.05 &
                                 abs(resOrdered_SN10_SS2$log2FoldChange)>=1,]
## This gave us 1804 genes that are differentially expressed using above criteria values

summary(res_clean_SN10_SS2)
## 21000 genes: 4323 genes upregulated and 4323 genes downregulated
sum(res_clean_SN10_SS2$padj < 0.1, na.rm=TRUE)
## 8646 genes differentially expressed
sum(res_clean_SN10_SS2$padj < 0.05, na.rm=TRUE)
## 7747 genes differentially expressed

S10_genes_ss2 <- as.character(sig_S10_SS2@rownames)


####################################################################################################################
## Differential expression analysis for technical replicates Sample S11 ##

setwd("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_2")
SN11_LPS_SS2 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_2/SN11_LPS_SS2.txt", quote="\"", comment.char="")
SN11_UNST_SS2 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_2/SN11_UNST_SS2.txt", quote="\"", comment.char="")

colnames (SN11_LPS_SS2) <- c("Gene_features", "LPS_Lane_4", "LPS_Lane_3", "LPS_Lane_2", "LPS_Lane_1")
colnames (SN11_UNST_SS2) <- c("Gene_features", "UNST_Lane_4", "UNST_Lane_3", "UNST_Lane_2", "UNST_Lane_1" ) 
data_SN11_SS2 <- merge (SN11_LPS_SS2, SN11_UNST_SS2, by.y ="Gene_features")
N = dim(data_SN11_SS2)[1]
rownames(data_SN11_SS2) = data_SN11_SS2[,1]
data_SN11_SS2 = data_SN11_SS2[,-1]
data_SN11_SS2 = data_SN11_SS2[c(6:N),]
data_SN11_SS2 <- data_SN11_SS2[ rowSums(data_SN11_SS2) > 1,]
## 26,209 genes in consideration

colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS", "LPS","UNST", "UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data_SN11_SS2, colData, formula(~ condition))
dds <- DESeq(dds)


pdf("~/BB2490-RNASeq-Project/results/MA_plots/MA_technical_replicates_SN11_SS2.pdf", height=8, width=12)
plotMA(dds, main="Differential Gene Expression in S11 technical replicates at Subsample 2")
dev.off()

## See the above argument to follow the steps done in the following steps
res_SN11_SS2 <- results(dds)
res_clean_SN11_SS2 <- res_SN11_SS2[(!is.na(res_SN11_SS2$padj)) &
                             (res_SN11_SS2$padj != 0.000000e+00), ]
resOrdered_SN11_ss2 <- res_clean_SN11_SS2[order(res_clean_SN11_SS2$padj),]
head(resOrdered_SN11_ss2)
sig_S11_SS2 <- resOrdered_SN11_ss2[resOrdered_SN11_ss2$padj<0.05 &
                                 abs(resOrdered_SN11_ss2$log2FoldChange)>=1,]
## This gave us 1043 genes that are differentially expressed using above criteria values

summary(res_clean_SN11_SS2)
## 3121 genes upregulated and 3519 genes downregulated
sum(res_clean_SN11_SS2$padj < 0.1, na.rm=TRUE)
## 6640 genes differentially expressed
sum(res_clean_SN11_SS2$padj < 0.05, na.rm=TRUE)
##  5827 genes differentially expressed

S11_genes_ss2 <- as.character(sig_S11_SS2@rownames)


####################################################################################################################
## Differential expression analysis for technical replicates Sample S12 ##

setwd("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_2")
SN12_LPS_SS2 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_2/SN12_LPS_SS2.txt", quote="\"", comment.char="")
SN12_UNST_SS2 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_2/SN12_UNST_SS2.txt", quote="\"", comment.char="")

colnames (SN12_LPS_SS2) <- c("Gene_features", "LPS_Lane_4", "LPS_Lane_3", "LPS_Lane_2", "LPS_Lane_1")
colnames (SN12_UNST_SS2) <- c("Gene_features", "UNST_Lane_4", "UNST_Lane_3", "UNST_Lane_2", "UNST_Lane_1") 
data_SN12_SS2 <- merge (SN12_LPS_SS2, SN12_UNST_SS2, by.y ="Gene_features")

N = dim(data_SN12_SS2)[1]
rownames(data_SN12_SS2) = data_SN12_SS2[,1]
data_SN12_SS2 = data_SN12_SS2[,-1]
data_SN12_SS2 = data_SN12_SS2[c(6:N),]
data_SN12_SS2 <- data_SN12_SS2[rowSums(data_SN12_SS2) > 1,]
## 26,536 genes in consideration

colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS", "LPS","UNST", "UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data_SN12_SS2, colData, formula(~ condition))
dds <- DESeq(dds)


pdf("~/BB2490-RNASeq-Project/results/MA_plots/MA_technical_replicates_SN12_SS2.pdf", height=8, width=12)
plotMA(dds, main="Differential Gene Expression in S12 technical replicates at Subsample 2")
dev.off()

## See the above argument to follow the steps done in the following steps
res_SN12_SS2 <- results(dds)
res_clean_SN12_SS2 <- res_SN12_SS2[(!is.na(res_SN12_SS2$padj)) &
                             (res_SN12_SS2$padj != 0.000000e+00), ]
resOrdered_SN12_SS2 <- res_clean_SN12_SS2[order(res_clean_SN12_SS2$padj),]
head(resOrdered_SN12_SS2)
sig_S12_SS2 <- resOrdered_SN12_SS2[resOrdered_SN12_SS2$padj<0.05 &
                                 abs(resOrdered_SN12_SS2$log2FoldChange)>=1,]
## This gave us 736 genes that are differentially expressed using above criteria values
## This gave us 833 genes that are differentially expressed using above criteria values

summary(resOrdered_SN12_SS2)
## 2107 genes upregulated and 2530 genes downregulated
## 2970 Upregulated and  3376  downregulated genes
sum(res_clean_SN12_SS2$padj < 0.1, na.rm=TRUE)
## 6346 genes differentially expressed
sum(res_clean_SN12_SS2$padj < 0.05, na.rm=TRUE)
##  5573 genes differentially expressed

S12_genes_ss2 <- as.character(sig_S12_SS2@rownames)


####################################################################################################################
####################################################################################################################
## Differential expression analysis for biological replicates all samples ##

Treated_Subsample_SS2 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_2/LPS_SS2.txt", quote="\"", comment.char="")
Untreated_Subsample_SS2 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_2/UNST_SS2.txt", quote="\"", comment.char="")

colnames (Untreated_Subsample_SS2) <- c("Gene_features", "UNST_SN10", "UNST_SN11", "UNST_SN12")
colnames (Treated_Subsample_SS2) <- c("Gene_features", "LPS_SN10", "LPS_SN11", "LPS_SN12") 
data_BS_SS2 <- merge (Treated_Subsample_SS2, Untreated_Subsample_SS2, by.y ="Gene_features" )
## Considering each lanes as technical replicates and hence given names based on lane numbers

N = dim(data_BS_SS2)[1]
rownames(data_BS_SS2) = data_BS_SS2[,1]
data_BS_SS2 = data_BS_SS2[,-1]
data_BS_SS2= data_BS_SS2[c(6:N),]
## removing last 5 rows which in our case turn out be in top 5 rows
data_BS_SS2 <- data_BS_SS2[ rowSums(data_BS_SS2) > 1, ] 
## Filtering to reduce number of genes that have 0 count values
## 27554 ENSEMBL genes
## 29784 ENSEMBL genes


colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS","UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data_BS_SS2, colData, formula(~ condition))
dds <- DESeq(dds)
# plotMA(dds, main="Differential Gene Expression in Sample S10 at Subsample 100% data")


pdf("~/BB2490-RNASeq-Project/results/MA_plots/MA_biological_replicates_SS2.pdf", height=8, width=12)
plotMA(dds, main="Differential Gene Expression in all samples at Subsample 2")
dev.off()
## The MA plot kinda of supports the argument of large number of genes differential expressed between two condition

res_SS2 <- results(dds)
res_clean_SS2 <- res_SS2[(!is.na(res_SS2$padj)) &
                   (res_SS2$padj != 0.000000e+00), ]
## I did this filtering to remove genes with 0 padjusted values
## May be it would be interesting to see why there is padjusted to 0

resOrdered_SS2 <- res_clean_SS2[order(res_clean_SS2$padj),]
head(resOrdered_SS2)
sig_SS2 <- resOrdered_SS2[resOrdered_SS2$padj<0.05 &
                    abs(resOrdered_SS2$log2FoldChange)>=1,]
## This gave us 1012 genes that are differential expressed with the above 
## criteria.
## Defning the criteria for the genes which are significant as ones
## with the Padjusted values lesser than 5% FDR and havin log2Fold change
## greater than 1. It would be interesting to see what happens with different
## cutoff. The above results ga

summary(res_clean_SS2)
#754 genes upregulated and 1247 genes downregulated
# 864 gene upregulated and 1449 genes downregulated
sum(res_clean_SS2$padj < 0.1, na.rm=TRUE)
## 2001 genes are differentially expressed at 10% FDR
## 2318 genes 

sum(res_clean_SS2$padj < 0.05, na.rm=TRUE)
## 1612 genes are differentially expressed at 5% FDR
## 1907 genes are differentially expressed at 5% FDR

genes_BR_SS2 <- as.character(sig_SS2@rownames)



#####################################
library(gplots)
library(VennDiagram)
Common_genes_SS2<- Reduce(intersect,  list(S10_genes_ss2, 
                                       S11_genes_ss2,
                                       S12_genes_ss2))

test_2 <- Reduce(intersect, list(Common_genes_SS2,
                               genes_BR_SS2))

#############################################################################################################

