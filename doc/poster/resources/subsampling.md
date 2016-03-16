We have developed a workflow which outputs a subsampling scheme based on a directory of BAM/SAM files and then generates subsamples with Samtools. Average read count is calculated for each sample group. Sample groups are then ordered by average read count. Our script does the following.

Step 1: A cut variable is defined in millions of reads. In Figure 1, 10 million was used.

Step 2: Whereever the difference between two average read counts is larger than the cut, the highest read count in the lowest of the two sample groups is defined as the maximum number of reads for that subsample. 

A smallest subsample cut can be defined by cutting so that no sample is fully included, if Figure 1 it is 50% of the reads in sample SN10 LPS L1.

———

When doing RNA-seq for multiple samples and technical replicates, sequencing depth can vary widely between samples. In this study, we try to investigate how subsampling of datasets unevenly distributed coverage affects differential expression analysis of  genes between untreated blood samples and samples treated with LPS to simulate inflammation. 

Dataset consists of Illunmina HiSeq paired end reads, 100 bp long. Blood samples from 3 individuals have been subjected to RNA-Seq. Two samples have been collected per individual, and one of each was treated with lipopolysaccharides (LPS) in order to stimulate an inflammatory response.


———

Differential expression analysis

Differential expression analysis was carried out comparing treated and untreated samples at FDR < 0.05% and at a fold change higher than 1 for all subsampled datasets. Genes meeting these requirements are thought to be differentially expressed.

Between subsamples 3 and 4 there is no significant change in differentially expressed genes. There is a small increase in DE genes between SS2 and SS3, and a substantial increase between SS1 and SS4.




