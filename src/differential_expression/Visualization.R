## Getting mapped reads from Whole genome and whole exome sequencing ###

install.packages("readxl")
library(readxl)
Read2features_SS1 <- read_excel("~/Desktop/Analysis of data from high throughput experiments/Features_Mapped_Count.xlsx",
                               sheet = "SS1")
Read2features_SS2 <-read_excel("~/Desktop/Analysis of data from high throughput experiments/Features_Mapped_Count.xlsx",
                               sheet = "SS2")
Read2features_SS3 <-read_excel("~/Desktop/Analysis of data from high throughput experiments/Features_Mapped_Count.xlsx",
                               sheet = "SS3")
Read2features_SS4 <-read_excel("~/Desktop/Analysis of data from high throughput experiments/Features_Mapped_Count.xlsx",
                               sheet = "SS4")


rownames(Read2features_SS1) = Read2features_SS1[,1]
rownames(Read2features_SS2) = Read2features_SS2[,1]
rownames(Read2features_SS3) = Read2features_SS3[,1]
rownames(Read2features_SS4) = Read2features_SS4[,1]

Read2features_SS1 = Read2features_SS1[,-1]
Read2features_SS2 = Read2features_SS2[,-1]
Read2features_SS3 = Read2features_SS3[,-1]
Read2features_SS4 = Read2features_SS4[,-1]


pdf("~/BB2490-RNASeq-Project/results/Count_Features.pdf", height=15, width=20)


Most
dev.off()


#############################################################################################################################

Differential_Expressed <- read_excel("~/BB2490-RNASeq-Project/results/Differential_Expression_Genes.xlsx",
                                sheet = "DE_Foldchange")


rownames(Differential_Expressed) = Differential_Expressed[,1]
Differential_Expressed = Differential_Expressed[,-1]


pdf("~/BB2490-RNASeq-Project/results/DE_genes.pdf", height=8, width=10)
barplot(as.matrix(Differential_Expressed), 
        ylim =c(0, max(Differential_Expressed)+120),
        beside=TRUE,
        col=c("dark green", "red", "blue"),
        ylab ="Counts of Differential Genes in each subsample",No
        main ="Differentually expressed genes at FDR < 0.05 and Foldchange >1 ",
        cex.axis=1, cex.lab=1.25, cex.main =1.25
)

legend("topleft",
       legend =  c("Upregulated Genes","Downregulated Genes","Total Differential Genes"), 
       col=c("dark green", "red", "blue"), lwd=5, horiz = TRUE,box.col = "white",bg = "white" )


dev.off()

##############################################################################################
library(gplots)
library(VennDiagram)
Common_genes<- Reduce(intersect,  list(sig@rownames, 
                                       sig_SS2@rownames,
                                       sig_SS3@rownames,
                                       sig_SS4@rownames))
## 824 genes 

venn.plot<-venn.diagram(x = list(SS1=sig@rownames, 
                                 SS2=sig_SS2@rownames,
                                 SS3=sig_SS3@rownames,
                                 SS4=sig_SS4@rownames),
                        filename = "Venn_2_plot.tiff",
                        col = "transparent",
                        fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
                        alpha = 0.50,
                        label.col = c("orange", "white", "darkorchid4", "white",
                                      "white", "white", "white", "white", "darkblue", "white",
                                      "white", "white", "white", "darkgreen", "white"),
                        cex = 1.5,
                        fontfamily = "serif",
                        fontface = "bold",
                        cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
                        cat.cex = 1.5,
                        cat.pos = 0,
                        cat.dist = 0.07,
                        cat.fontfamily = "serif",
                        rotation.degree = 0,
                        main =" Venn Diagram depicting common DE identified in each subsample"
)

