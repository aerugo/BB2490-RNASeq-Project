##Day -1##
Log book

##Day 2:##

#Project Title and intial feel of the project#

First Observation: Title : The effect of sub-sampling of large RNA-seq data sets on gene expression analysis

Keywords: Subsampling, large RNAseq, gene expression analysis

#Introduction:#

Diagnosis the problem:
*Normal blood cells from individuals were collected and a part treated with LPS while other untreated
*A high depth sequencing experiments performed in these LPS treated and normal sample
*Rather few gene expression observed. Hypothetically, a high sequencing leading to inability to identify differentional status 

#Expectation#
*Differential expressed and Allele specific genes should be repored using approproate methods and guidelines

#Methods and Materials#
*Subsampling of Reads
*30 - 70% of data take into consideration and perform the DE and ASE. 
*2 folds improvements --> less reads to be aligned and more more meaningful interpretstion?

What is subsampling and how can you make sampling such that it represents the true population?
Selection of the sample such that it characterizes population.

How are previous study done in in this case? 
Why does sampling matters ?
How do you proceed?
take 30%,40%,50%,60%,70% of reads and do Trimming, QC, Alignment, DE, ASE of the reads?
 8 samples X 5 times run the pipelines?
 
Pipleines:
*Preprocessing : 
**Step 1: Removal of adapters, low quality reads with q values less than some parameters; 
*Tools in consideration -->Trimgalore, Cutadapt
Step 2: Quality Control of the trimmed reads 
Tools in consideration FASTQC

Sanity check

#Alignment and Mapping#

Various mappers BWA, bowtie2, STAR, Tophat, But given the RNASeq data we have to consider take into consideration of spliced reads and computation time for the mapping the reads

QC of the mapping 
How good are our read maps?
Where do they align?
How much of the genome is covered? 
Spliced site informations?

Tools: HS metrices, Samtools 

Differential gene expression analysis:
Find the expression fold change and p value of the differential expressed genes
Does this make sense?
How would the subsampling effect the differential expressed genes efect?

Allele specific expression 
What is allele specific expression?
Allele specific expression is the difference in the paternal and maternal haplotype from an individual. ## Very very intersting  
