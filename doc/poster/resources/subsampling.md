We have developed a workflow which outputs a subsampling scheme based on a directory of BAM/SAM files and then generates subsamples with Samtools. Average read count is calculated for each sample group. Sample groups are then ordered by average read count. Our script does the following.

Step 1: A cut variable is defined in millions of reads. In Figure 1, 10 million was used.

Step 2: Whereever the difference between two average read counts is larger than the cut, the highest read count in the lowest of the two sample groups is defined as the maximum number of reads for that subsample. 

A smallest subsample cut can be defined by cutting so that no sample is fully included, if Figure 1 it is 50% of the reads in sample SN10 LPS L1.

