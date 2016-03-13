## Finally while doing the subsampling scheme from SS1 to SS4 only LPS SS1 have enough of reads to contribute to the subsampling 4 . Hence only
## differential expression analysis of Sample S11 is done .

## join -j 1 SN_11_LPS_TGACCA_L003_count.txt SN_11_LPS_TGACCA_L002_count.txt| join -j 1 SN_11_LPS_TGACCA_L004_count.txt -  | join -j 1 SN_11_LPS_TGACCA_L005_count.txt - > SN11_LPS_SS4.txt


## awk '{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print $1 " " sum}' SN11_LPS_SS4.txt > LPS_temp2

# join -j 1  LPS_temp2 LPS_temp3| join -j 1 LPS_temp1 -  > LPS_SS3.txt


############################################################################################
## Rest of the data are imputed from Subsampling scheme SS2 and SS3 respectively for other samples

## Differential expression analysis for technical replicates Sample S11 ##
setwd("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_4")
SN11_LPS_SS4 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_4/SN11_LPS_SS4.txt", quote="\"", comment.char="")
SN11_UNST_SS4 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_4/SN11_UNST_SS4.txt", quote="\"", comment.char="")

colnames (SN11_LPS_SS4) <- c("Gene_features", "LPS_Lane_4", "LPS_Lane_3", "LPS_Lane_2", "LPS_Lane_1")
colnames (SN11_UNST_SS4) <- c("Gene_features", "UNST_Lane_4", "UNST_Lane_3", "UNST_Lane_2", "UNST_Lane_1" ) 
data_SN11_SS4 <- merge (SN11_LPS_SS4, SN11_UNST_SS4, by.y ="Gene_features")
N = dim(data_SN11_SS4)[1]
rownames(data_SN11_SS4) = data_SN11_SS4[,1]
data_SN11_SS4 = data_SN11_SS4[,-1]
data_SN11_SS4 = data_SN11_SS4[c(6:N),]
data_SN11_SS4 <- data_SN11_SS4[ rowSums(data_SN11_SS4) > 1,]
## 29,307 genes in consideration

colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS", "LPS","UNST", "UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data_SN11_SS4, colData, formula(~ condition))
dds <- DESeq(dds)


pdf("~/BB2490-RNASeq-Project/results/MA_plots/MA_technical_replicates_SN11_SS4.pdf", height=8, width=12)
plotMA(dds, main="Differential Gene Expression in S11 techincal replicates at Subsample 4")
dev.off()

## See the above argument to follow the steps done in the following steps
res_SN11_SS4 <- results(dds)
res_clean_SN11_SS4 <- res_SN11_SS4[(!is.na(res_SN11_SS4$padj)) &
                             (res_SN11_SS4$padj != 0.000000e+00), ]
resOrdered_SN11_SS4 <- res_clean_SN11_SS4[order(res_clean_SN11_SS4$padj),]
head(resOrdered_SN11_SS4)
sig_S11_SS4 <- resOrdered_SN11_SS4[resOrdered_SN11_SS4$padj<0.05 &
                                 abs(resOrdered_SN11_SS4$log2FoldChange)>=1,]
## This gave us 1169 genes that are differentially expressed using above criteria values

summary(res_clean_SN11_SS4)
## 4068 genes upregulated and 4367 genes downregulated
sum(res_clean_SN11_SS4$padj < 0.1, na.rm=TRUE)
## 8435) genes differentially expressed
sum(res_clean_SN11_SS4$padj < 0.05, na.rm=TRUE)
##  7557 genes differentially expressed

S11_genes_ss4 <- as.character(sig_S11_SS4@rownames)



####################################################################################################################
####################################################################################################################
## Differential expression analysis for biological replicates all samples ##

Treated_Subsample_SS4 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_4/LPS_SS4.txt", quote="\"", comment.char="")
Untreated_Subsample_SS4 <- read.table("~/BB2490-RNASeq-Project/data/count/Htseq/Subsample_4/UNST_SS4.txt", quote="\"", comment.char="")

colnames (Untreated_Subsample_SS4) <- c("Gene_features", "UNST_SN10", "UNST_SN11", "UNST_SN12")
colnames (Treated_Subsample_SS4) <- c("Gene_features", "LPS_SN10", "LPS_SN11", "LPS_SN12") 
data_BS_SS4 <- merge (Treated_Subsample_SS4, Untreated_Subsample_SS4, by.y ="Gene_features" )
## Considering each lanes as technical replicates and hence given names based on lane numbers

N = dim(data_BS_SS4)[1]
rownames(data_BS_SS4) = data_BS_SS4[,1]
data_BS_SS4 = data_BS_SS4[,-1]
data_BS_SS4 = data_BS_SS4[c(6:N),]
## removing last 5 rows which in our case turn out be in top 5 rows
data_BS_SS4 <- data_BS_SS4[ rowSums(data_BS_SS4) > 1, ] 
## Filtering to reduce number of genes that have 0 count values
## 27554 ENSEMBL genes
## 29784 ENSEMBL genes
## 31106 ENSEMBL genes
## 31594 ENSEMBL genes


colData <- DataFrame(condition=factor(c("LPS", "LPS", "LPS","UNST", "UNST", "UNST")))
dds <- DESeqDataSetFromMatrix(data_BS_SS4, colData, formula(~ condition))
dds <- DESeq(dds)
# plotMA(dds, main="Differential Gene Expression in Sample S10 at Subsample 100% data")


pdf("~/BB2490-RNASeq-Project/results/MA_plots/MA_biological_replicates_SS4.pdf", height=8, width=12)
plotMA(dds, main="Differential Gene Expression in all samples at Subsample 4")
dev.off()
## The MA plot kinda of supports the argument of large number of genes differential expressed between two condition

res_SS4 <- results(dds)
res_clean_SS4 <- res_SS4[(!is.na(res_SS4$padj)) &
                           (res_SS4$padj != 0.000000e+00), ]
## I did this filtering to remove genes with 0 padjusted values
## May be it would be interesting to see why there is padjusted to 0

resOrdered_SS4 <- res_clean_SS4[order(res_clean_SS4$padj),]
head(resOrdered_SS4)
sig_SS4 <- resOrdered_SS4[resOrdered_SS4$padj<0.05 &
                            abs(resOrdered_SS4$log2FoldChange)>=1,]
## This gave us 1076 genes that are differential expressed with the above 
## criteria.
## Defning the criteria for the genes which are significant as ones
## with the Padjusted values lesser than 5% FDR and havin log2Fold change
## greater than 1. It would be interesting to see what happens with different
## cutoff. The above results ga

summary(res_clean_SS4)
#754 genes upregulated and 1247 genes downregulated
# 864 gene upregulated and 1449 genes downregulated
# 954 gene upregulated and 1569 downregulated
# 980 genes upregulated aand 1600 downregulated
sum(res_clean_SS4$padj < 0.1, na.rm=TRUE)
## 2001 genes are differentially expressed at 10% FDR
## 2318 genes 
##  2527
## 2580

sum(res_clean_SS4$padj < 0.05, na.rm=TRUE)
## 1612 genes are differentially expressed at 5% FDR
## 1907 genes are differentially expressed at 5% FDR
## 2053
## 2109

genes_BR_SS4 <- as.character(sig_SS4@rownames)
