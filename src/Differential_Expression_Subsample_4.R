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

## Differential expression analysis for Sample S10 ##

setwd("~/BB2490-RNASeq-Project/data/Htseq")
SN10_LPS_Subsample4 <- read.table("~/BB2490-RNASeq-Project/data/Htseq/SN10_LPS_Subsample4.txt", quote="\"", comment.char="")
SN10_UNST_Subsample4 <- read.table("~/BB2490-RNASeq-Project/data/Htseq/SN10_UNST_Subsample4.txt", quote="\"", comment.char="")

colnames (SN10_UNST_Subsample4) <- c("Gene_features", "UNST_Lane_1", "UNST_Lane_2", "UNST_Lane_3", "UNST_Lane_4")
colnames (SN10_LPS_Subsample4) <- c("Gene_features", "LPS_Lane_1", "LPS_Lane_2", "LPS_Lane_3", "LPS_Lane_4") 
data <- merge (SN10_LPS_Subsample4, SN10_UNST_Subsample4, by.y ="Gene_features"  )

N = dim(data)[1]
rownames(data) = data[,1]
data = data[,-1]
data = data[c(6:N),]
data <- data[ rowSums(data) > 1, ]

colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS", "LPS","UNST", "UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data, colData, formula(~ condition))
dds <- DESeq(dds)
plotMA(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
head(resOrdered)
sig <- resOrdered[!is.na(resOrdered$padj) &
                    resOrdered$padj<0.05 &
                    abs(resOrdered$log2FoldChange)>=1,]
True_diff<-sig[which(sig$padj != 0.000000e+00),]

summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
sum(res$padj < 0.01, na.rm=TRUE)



## Differential expresssion concencerning 100% of data from flowcells

## DEQSeq2 : Differential gene expression analysis in Inflammatoy and Non Inflammatory agents in WBC ##
library("DESeq2")
## Differential expression analysis for Sample S11 ##

setwd("~/BB2490-RNASeq-Project/data/Htseq")
SN11_LPS_Subsample4 <- read.table("~/BB2490-RNASeq-Project/data/Htseq/SN11_LPS_Subsample4.txt", quote="\"", comment.char="")
SN11_UNST_Subsample4 <- read.table("~/BB2490-RNASeq-Project/data/Htseq/SN11_UNST_Subsample4.txt", quote="\"", comment.char="")

colnames (SN11_LPS_Subsample4) <- c("Gene_features", "LPS_Lane_1", "LPS_Lane_2", "LPS_Lane_3", "LPS_Lane_4")
colnames (SN11_UNST_Subsample4) <- c("Gene_features", "UNST_Lane_1", "UNST_Lane_2", "UNST_Lane_3", "UNST_Lane_4" ) 
data <- merge (SN11_LPS_Subsample4, SN11_UNST_Subsample4, by.y ="Gene_features"  )

N = dim(data)[1]
rownames(data) = data[,1]
data = data[,-1]
data = data[c(6:N),]
data <- data[ rowSums(data) > 1, ]

colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS", "LPS","UNST", "UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data, colData, formula(~ condition))
dds <- DESeq(dds)
plotMA(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
head(resOrdered)
sig <- resOrdered[!is.na(resOrdered$padj) &
                    resOrdered$padj<0.05 &
                    abs(resOrdered$log2FoldChange)>=1,]
True_diff<-sig[which(sig$padj != 0.000000e+00),]


summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
sum(res$padj < 0.01, na.rm=TRUE)


## Differential expression analysis for Sample S12 ##

setwd("~/BB2490-RNASeq-Project/data/Htseq")
SN12_LPS_Subsample4 <- read.table("~/BB2490-RNASeq-Project/data/Htseq/SN12_LPS_Subsample4.txt", quote="\"", comment.char="")
SN12_UNST_Subsample4 <- read.table("~/BB2490-RNASeq-Project/data/Htseq/SN12_UNST_Subsample4.txt", quote="\"", comment.char="")

colnames (SN12_LPS_Subsample4) <- c("Gene_features", "LPS_Lane_1", "LPS_Lane_2", "LPS_Lane_3", "LPS_Lane_4")
colnames (SN12_UNST_Subsample4) <- c("Gene_features", "UNST_Lane_1", "UNST_Lane_2", "UNST_Lane_3", "UNST_Lane_4" ) 
data <- merge (SN12_LPS_Subsample4, SN12_UNST_Subsample4, by.y ="Gene_features"  )

N = dim(data)[1]
rownames(data) = data[,1]
data = data[,-1]
data = data[c(6:N),]
data <- data[ rowSums(data) > 1, ]

colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS", "LPS","UNST", "UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data, colData, formula(~ condition))
dds <- DESeq(dds)
plotMA(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
head(resOrdered)
sig <- resOrdered[!is.na(resOrdered$padj) &
                    resOrdered$padj<0.05 &
                    abs(resOrdered$log2FoldChange)>=1,]
True_diff<-sig[which(sig$padj != 0.000000e+00),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
sum(res$padj < 0.01, na.rm=TRUE)

