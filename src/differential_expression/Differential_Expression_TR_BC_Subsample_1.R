## Differential Expression analysis for Subsample 1 

###### Awk codes for combining files together  for technical replicates (Possible for cooler solutions)
## join -j 1  SN_10_LPS_CGATGT_L002_count.txt SN_10_LPS_CGATGT_L001_count.txt| join -j 1 SN_10_LPS_CGATGT_L003_count.txt -  | join -j 1 SN_10_LPS_CGATGT_L004_count.txt - > SN10_LPS_SS1.txt
## join -j 1  SN10_UNST_ATCACG_L002_count.txt SN10_UNST_ATCACG_L001_count.txt| join -j 1 SN10_UNST_ATCACG_L003_count.txt -  | join -j 1 SN10_UNST_ATCACG_L004_count.txt - > SN10_UNST_SS1.txt
## join -j 1  SN11_UNST_TTAGGC_L003_count.txt SN11_UNST_TTAGGC_L002_count.txt| join -j 1 SN11_UNST_TTAGGC_L004_count.txt -  | join -j 1 SN11_UNST_TTAGGC_L005_count.txt - > SN11_UNST_SS1.txt
## join -j 1  SN_11_LPS_TGACCA_L003_count.txt SN_11_LPS_TGACCA_L002_count.txt| join -j 1 SN_11_LPS_TGACCA_L004_count.txt -  | join -j 1 SN_11_LPS_TGACCA_L005_count.txt - > SN11_LPS_SS1.txt
## join -j 1  SN_12_LPS_GCCAAT_L004_count.txt SN_12_LPS_GCCAAT_L003_count.txt| join -j 1 SN_12_LPS_GCCAAT_L005_count.txt -  | join -j 1 SN_12_LPS_GCCAAT_L006_count.txt - > SN12_LPS_SS1.txt
## join -j 1  SN12_UNST_ACAGTG_L004_count.txt SN12_UNST_ACAGTG_L003_count.txt| join -j 1 SN12_UNST_ACAGTG_L005_count.txt -  | join -j 1 SN12_UNST_ACAGTG_L006_count.txt - > SN12_UNST_SS1.txt

## For Biological Replicates
# awk '{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print $1 " " sum}' SN10_LPS_SS1.txt > LPS_temp1
# awk '{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print $1 " " sum}' SN11_LPS_SS1.txt > LPS_temp2
# awk '{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print $1 " " sum}' SN12_LPS_SS1.txt > LPS_temp3
# 
# awk '{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print $1 " " sum}' SN10_UNST_SS1.txt > UNST_temp1
# awk '{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print $1 " " sum}' SN11_UNST_SS1.txt > UNST_temp2
# awk '{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print $1 " " sum}' SN12_UNST_SS1.txt > UNST_temp3
# 
# join -j 1  LPS_temp2 LPS_temp3| join -j 1 LPS_temp1 -  > LPS_SS1.txt
# join -j 1  UNST_temp2 UNST_temp3| join -j 1 UNST_temp1 -  > UNST_SS1.txt

####################################################################################################################
## Differential expression for technical replicates considering 

library("DESeq2")
####################################################################################################################
## Differential expression analysis for Sample S10 ##

setwd("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_1")
SN10_LPS_SS1 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_1/SN10_LPS_SS1.txt", quote="\"", comment.char="")
SN10_UNST_SS1 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_1/SN10_UNST_SS1.txt", quote="\"", comment.char="")

colnames (SN10_LPS_SS1) <- c("Gene_features", "LPS_Lane_4", "LPS_Lane_3", "LPS_Lane_2", "LPS_Lane_1")
colnames (SN10_UNST_SS1) <- c("Gene_features", "UNST_Lane_4", "UNST_Lane_3", "UNST_Lane_2", "UNST_Lane_1" ) 
data_SN10_SS1 <- merge (SN10_LPS_SS1, SN10_UNST_SS1, by.y ="Gene_features")

N = dim(data_SN10_SS1)[1]
rownames(data_SN10_SS1) = data_SN10_SS1[,1]
data_SN10_SS1 = data_SN10_SS1[,-1]
data_SN10_SS1 = data_SN10_SS1[c(6:N),]
data_SN10_SS1 <- data_SN10_SS1[ rowSums(data_SN10_SS1) > 1,]
## 23,711 genes in consideration

colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS", "LPS","UNST", "UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data_SN10_SS1, colData, formula(~ condition))
dds <- DESeq(dds)


pdf("~/BB2490-RNASeq-Project/results/MA_plots/MA_technical_replicates_SN10_SS1.pdf", height=8, width=12)
plotMA(dds, main="Differential Gene Expression in S10 techincal replicates at Subsample 1")
dev.off()

## See the above argument to follow the steps done in the following steps
res_SN10 <- results(dds)
res_clean_SN10 <- res_SN10[(!is.na(res_SN10$padj)) &
                             (res_SN10$padj != 0.000000e+00), ]
resOrdered_SN10 <- res_clean_SN10[order(res_clean_SN10$padj),]
head(resOrdered_SN10)
sig_S10_SS1 <- resOrdered_SN10[resOrdered_SN10$padj<0.05 &
                             abs(resOrdered_SN10$log2FoldChange)>=1,]
## This gave us 1671 genes that are differentially expressed using above criteria values

summary(res_clean_SN10)
## 3559 genes upregulated and 3613 genes downregulated
sum(res_clean_SN10$padj < 0.1, na.rm=TRUE)
## 7172 genes differentially expressed
sum(res_clean_SN10$padj < 0.05, na.rm=TRUE)
## 6386 genes differentially expressed

S10_genes_ss1 <- as.character(sig_S10_SS1@rownames)


####################################################################################################################
## Differential expression analysis for technical replicates Sample S11 ##

setwd("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_1")
SN11_LPS_SS1 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_1/SN11_LPS_SS1.txt", quote="\"", comment.char="")
SN11_UNST_SS1 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_1/SN11_UNST_SS1.txt", quote="\"", comment.char="")

colnames (SN11_LPS_SS1) <- c("Gene_features", "LPS_Lane_4", "LPS_Lane_3", "LPS_Lane_2", "LPS_Lane_1")
colnames (SN11_UNST_SS1) <- c("Gene_features", "UNST_Lane_4", "UNST_Lane_3", "UNST_Lane_2", "UNST_Lane_1" ) 
data_SN11_SS1 <- merge (SN11_LPS_SS1, SN11_UNST_SS1, by.y ="Gene_features")

N = dim(data_SN11_SS1)[1]
rownames(data_SN11_SS1) = data_SN11_SS1[,1]
data_SN11_SS1 = data_SN11_SS1[,-1]
data_SN11_SS1 = data_SN11_SS1[c(6:N),]
data_SN11_SS1 <- data_SN11_SS1[ rowSums(data_SN11_SS1) > 1,]
## 23,824 genes in consideration

colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS", "LPS","UNST", "UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data_SN11_SS1, colData, formula(~ condition))
dds <- DESeq(dds)


pdf("~/BB2490-RNASeq-Project/results/MA_plots/MA_technical_replicates_SN11_SS1.pdf", height=8, width=12)
plotMA(dds, main="Differential Gene Expression in S11 techincal replicates at Subsample 1")
dev.off()

## See the above argument to follow the steps done in the following steps
res_SN11 <- results(dds)
res_clean_SN11 <- res_SN11[(!is.na(res_SN11$padj)) &
                             (res_SN11$padj != 0.000000e+00), ]
resOrdered_SN11 <- res_clean_SN11[order(res_clean_SN11$padj),]
head(resOrdered_SN11)
sig_S11_SS1 <- resOrdered_SN11[resOrdered_SN11$padj<0.05 &
                                 abs(resOrdered_SN11$log2FoldChange)>=1,]
## This gave us 917 genes that are differentially expressed using above criteria values

summary(res_clean_SN11)
## 2302 genes upregulated and 2703 genes downregulated
sum(res_clean_SN11$padj < 0.1, na.rm=TRUE)
## 5005 genes differentially expressed
sum(res_clean_SN11$padj < 0.05, na.rm=TRUE)
##  4374 genes differentially expressed

S11_genes_ss1 <- as.character(sig_S11_SS1@rownames)


####################################################################################################################
## Differential expression analysis for technical replicates Sample S12 ##

setwd("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_1")
SN12_LPS_SS1 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_1/SN12_LPS_SS1.txt", quote="\"", comment.char="")
SN12_UNST_SS1 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_1/SN12_UNST_SS1.txt", quote="\"", comment.char="")

colnames (SN12_LPS_SS1) <- c("Gene_features", "LPS_Lane_4", "LPS_Lane_3", "LPS_Lane_2", "LPS_Lane_1")
colnames (SN12_UNST_SS1) <- c("Gene_features", "UNST_Lane_4", "UNST_Lane_3", "UNST_Lane_2", "UNST_Lane_1" ) 
data_SN12_SS1 <- merge (SN12_LPS_SS1, SN12_UNST_SS1, by.y ="Gene_features")

N = dim(data_SN12_SS1)[1]
rownames(data_SN12_SS1) = data_SN12_SS1[,1]
data_SN12_SS1 = data_SN12_SS1[,-1]
data_SN12_SS1 = data_SN12_SS1[c(6:N),]
data_SN12_SS1 <- data_SN12_SS1[rowSums(data_SN12_SS1) > 1,]
## 23,987 genes in consideration

colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS", "LPS","UNST", "UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data_SN12_SS1, colData, formula(~ condition))
dds <- DESeq(dds)


pdf("~/BB2490-RNASeq-Project/results/MA_plots/MA_technical_replicates_SN12_SS1.pdf", height=8, width=12)
plotMA(dds, main="Differential Gene Expression in S12 techincal replicates at Subsample 1")
dev.off()

## See the above argument to follow the steps done in the following steps
res_SN12 <- results(dds)
res_clean_SN12 <- res_SN12[(!is.na(res_SN12$padj)) &
                             (res_SN12$padj != 0.000000e+00), ]
resOrdered_SN12 <- res_clean_SN12[order(res_clean_SN12$padj),]
head(resOrdered_SN12)
sig_S12_SS1 <- resOrdered_SN12[resOrdered_SN12$padj<0.05 &
                                 abs(resOrdered_SN12$log2FoldChange)>=1,]
## This gave us 736 genes that are differentially expressed using above criteria values

summary(res_clean_SN12)
## 2107 genes upregulated and 2530 genes downregulated
sum(res_clean_SN12$padj < 0.1, na.rm=TRUE)
## 4637 genes differentially expressed
sum(res_clean_SN12$padj < 0.05, na.rm=TRUE)
##  3981 genes differentially expressed

S12_genes_ss1 <- as.character(sig_S12_SS1@rownames)


####################################################################################################################
####################################################################################################################
## Differential expression analysis for biological replicates Sample S12 ##

Treated_Subsample <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_1/LPS_SS1.txt", quote="\"", comment.char="")
Untreated_Subsample <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_1/UNST_SS1.txt", quote="\"", comment.char="")

colnames (Untreated_Subsample) <- c("Gene_features", "UNST_SN10", "UNST_SN11", "UNST_SN12")
colnames (Treated_Subsample) <- c("Gene_features", "LPS_SN10", "LPS_SN11", "LPS_SN12") 
data_BS <- merge (Treated_Subsample, Untreated_Subsample, by.y ="Gene_features" )
## Considering each lanes as technical replicates and hence given names based on lane numbers

N = dim(data_BS)[1]
rownames(data_BS) = data_BS[,1]
data_BS = data_BS[,-1]
data_BS= data_BS[c(6:N),]
## removing last 5 rows which in our case turn out be in top 5 rows
data_BS <- data_BS[ rowSums(data_BS) > 1, ] 
## Filtering to reduce number of genes that have 0 count values
## 27554 ENSEMBL genes

colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS","UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data_BS, colData, formula(~ condition))
dds <- DESeq(dds)
# plotMA(dds, main="Differential Gene Expression in Sample S10 at Subsample 100% data")


pdf("~/BB2490-RNASeq-Project/results/MA_plots/MA_biological_replicates_SS1.pdf", height=8, width=12)
plotMA(dds, main="Differential Gene Expression in all samples at Subsample 1")
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
## This gave us 870 genes that are differential expressed with the above 
## criteria.
## Defning the criteria for the genes which are significant as ones
## with the Padjusted values lesser than 5% FDR and havin log2Fold change
## greater than 1. It would be interesting to see what happens with different
## cutoff. The above results ga

summary(res_clean)
#754 genes upregulated and 1247 genes downregulated
sum(res_clean$padj < 0.1, na.rm=TRUE)
## 2001 genes are differentially expressed at 10% FDR

sum(res_clean$padj < 0.05, na.rm=TRUE)
## 1612 genes are differentially expressed at 1% FDR

genes_BR <- as.character(sig@rownames)



#####################################
library(gplots)
library(VennDiagram)
Common_genes<- Reduce(intersect,  list(S10_genes_ss1, 
                                       S11_genes_ss1,
                                       S12_genes_ss1))

test <- Reduce(intersect, list(Common_genes,
                               genes_BR))

#############################################################################################################

