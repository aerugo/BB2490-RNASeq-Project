#Day -1#
Log book

#Day 2: Feb 18 2016 ##

##Project Title and intial feel of the project##

First Observation: Title : The effect of sub-sampling of large RNA-seq data sets on gene expression analysis

Keywords: Subsampling, large RNAseq, gene expression analysis

###Introduction:###

Diagnosis the problem:
* Normal blood cells from individuals were collected and a part treated with LPS while other untreated
* A high depth sequencing experiments performed in these LPS treated and normal sample
* Rather few gene expression observed. Hypothetically, a high sequencing leading to inability to identify differentional status 

###Aim####
* Differential expressed and Allele specific genes should be repored using approproate methods and guidelines

###Methods and Materials####
* Subsampling of Reads
* 30 - 70% of data take into consideration and perform the DE and ASE. 
* 2 folds improvements --> less reads to be aligned and more more meaningful interpretstion?

What is subsampling and how can you make sampling such that it represents the true population?
* Selection of the sample such that it characterizes population.

* How are previous study done in in this case? 
* Why does sampling matters ?
* How do you proceed?
 take 30%,40%,50%,60%,70% of reads and do Trimming, QC, Alignment, DE, ASE of the reads?
* 8 samples X 5 times run the pipelines?
 
Pipleines:
* Preprocessing : 
* Step 1: Removal of adapters, low quality reads with q values less than some parameters; 
* Tools in consideration -->Trimgalore, Cutadapt
* Step 2: Quality Control of the trimmed reads 
* Tools in consideration FASTQC

Sanity check

###Alignment and Mapping###

* Various mappers BWA, bowtie2, STAR, Tophat, 
* But given the RNASeq data we have to consider take into consideration of spliced reads and computation time for the mapping the reads

* QC of the mapping 
* How good are our read maps?
* Where do they align?
* How much of the genome is covered? 
* Spliced site informations?

Tools: HS metrices, Samtools 

* Differential gene expression analysis:
* Find the expression fold change and p value of the differential expressed genes
* Does this make sense?
* How would the subsampling effect the differential expressed genes efect?

Allele specific expression 
* What is allele specific expression?
* Allele specific expression is the difference in the paternal and maternal haplotype from an individual. ## Very very intersting  

# Day- 3: Feb 19 2016#
We further more discussed the project in detail with teacher regarding basic ideas and motive of project. Intailly, I felt the project is vague in sense that my major concern is as we sequence more the better would be results be. However, due to increased sequencing, it kinda of obsucures the major signal. Thus a subsampling technique, which would in principle implicate similar properties as population will provide better signal to noise ratio. This would in turn lead to better identification of gene expression and allele specific expression. 

Furthermore, I was confused about allele specific expressions. Through discussion in my own words, allele specific expression can be described as expression of paternal and maternal allele in an inddividual. Basically, an individual would have two set of chromosome; each arising from father and mother and in genral case a haloptype would be build based on the above information. However, I have to figure out how could paternal and maternal side of expression would be determined from RNA-seq data.

## General experimental design from RNAseq ##:
* 3 sample of each such that there are 6 replicates
* 2 group for Differential expression
* 6 individual sample with allele specific expression.

# Day 3 and 4 : 20/21 Feb 2016##

##Literature review and brief summary of each papers##

[Effects of subsampling on characteristics of RNA-seq data from triple-negative breast cancer patients][1] (http://cjcjournal.biomedcentral.com/articles/10.1186/s40880-015-0040-8)

* Keypoints :
* Subsampling of RNA reads represented a biological simulation albet in low coverage depth
* Saturation result -> 32 Million average sequencing depth and 46 million reads for an individual
* Moment in statistics would mean quantitative measures. Higher moment of expressed genes are better measures than mean of expressed genes

Triple negative breast cancer cells (TNBC) and triple negative breast cancer free (TNBC-free) samples RNA sequenced data were analyzed. Each samples' read were repeatedly subsampled to on mapped read to obtained subsamples with different sequencing depths.RNA-seq count distibution for expressed genes were studied. A subsampling technique called as Depth of Sequencing Iterative Reduction Estimator
(DESIRE) was introduced.

###Basic principle of DESIRE ###
Systematically, bootstapping using *f* fraction of read data from whole dataset. Here, the data are drawn without replacement.For each sample a fraction ranging from 10%-90% of reads are drawn without replacement to give R replicates. R=24 taken on the basis of number of nodes in cluster. However, it seems it is kinda of optimal value.

[GeneiASE: Detection of condition-dependent and static allele-specific expression from RNA-seq data without haplotype information][2] (http://www.nature.com/articles/srep21134)

Talked about GeneiAse method for the detection of allele specific expression from RNAseq data from LPS treated and untreated normal blood cells. 

# Day 5: Getting hands on data #

As the teacher had set data, I first wanted to get inhand with RNA seq data. My aim was to understand the reads count in each of sample. Basically, the data consists of three RNA sequenced samples- The samples were numbered as 10, 11 and 12 and consists of normal lympohodial cell treated with LPS and without LPS. Thus consisting of six samples in total. Furthermore, the data consists of each 2 flowcells thus consisting of 8 lanes: 4 lanes of each treated and 4 lanes without LPS treated. My first task was to observe the read counts in each sample. For this, I first made a noback files in proj g2015056. 

```shell
mkdir /proj/g2015056/nobackup/Subsampling
```
And using unix command as below count the individiual reads in each fastq files

```shell
 zgrep -c '^@' Sample_SN_11_UNST/S* > /proj/g2015056/nobackup/Subsample/130104_SN866_0198_BC1DAYACXX/RC_SN11_UNST
```

Here we see quite a drastic range of sequence reads in each of files. For example in flowcell "130104_SN866_0197_AC1DLVACXX" Sample10 LPS range 1 was ranged from 16.95M to 21.20M range while in other flowcell it ranged from 17.38M to 21.06M. We plotted a barplot of reads from 2 flow cells as shown in figure (dont know how to post figure in here)

# Day 6: Working together and hands on sbatch script with Hugio#

As planned I and Hugio met at KTH to get hands on data analysis and substantiate pipelines for further analysis. Firstly, we decided to used step-wise increment of data. Specifically we decided first we will try with a single lane reads so as to optimize tools and sbatch script and subsequently scale it up with serially with two, three and four lanes of each flowcells thus constituting 25%, 50%, 75% and 100% of total reads. The main rationale behind it was this would be in principle similar to subsampling whole population of paired ends reads from whole RNA sequencing albeit in small chunks.

## Developing pipelines for analysis ##

In my opinion major task in sequencing is composed of three major parts :
* Preprocessing of reads
* Alignment and Mapping of reads
* Downstreaming analysis from mapped reads

Thus, in each step we have to rigorously determine tool's advantages and shortcomings . At first we started tentatively with listing (for preprocessing used) tools for analysis of RNA seq data.

* Preprocessing of read data

Tools in consideration : [TrimGalore] (http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_User_Guide_v0.3.7.pdf), [FastQC] (https://wiki.hpcc.msu.edu/display/Bioinfo/FastQC+Tutorial), [Trimmomatic] (http://bioinformatics.oxfordjournals.org/content/30/15/2114.full.pdf+html), [FastX] (http://hannonlab.cshl.edu/fastx_toolkit/index.html)













