# Log for BB2490 Project

## February 17th 2016
We settled on GitHub both for hosting code and logbooks. Since the rest of the group were not very familiar with Git, we had a little workshop this afternoon where we set up the GitHub repo. We also decided on Python as our go-to scripting language and set up the IDE PyCharm on Amy's and Tobias's computers. We also made sure that UppMax has Git installed. 

## February 19th 2016
Meeting with Olof after seminar. 
To reduce memory footprint and computational time we will restrict to 3 individuals, rather than 8.
First step will be to remove adapters, check quality and other preprocessing. 
We discussied how we would do subsampling. Current idea is that sampling will be done after processing but before mapping.
Each sample was sequenced on 8 different lanes, four from each flow cell.

## February 23d 2016
I worked on preprocessing with Sailendra this afternoon. We started by becoming familiar [with the data](https://github.com/aerugo/BB2490-RNASeq-Project/wiki/Overview-of-data). It is structured as follows:

```
|
|-- Flowcell 1: 130104_SN866_0197_AC1DLVACXX
|   |-- SN_10_LPS: Sample from person 1, treated with LPS
|   |   `-- Paired end reads in fastq format from four lanes of flowcell 1, adapter with tag CGATGT.
|   |-- SN_10_UNST: Sample from person 1, untreated
|   |   `-- Paired end reads in fastq format from four lanes of flowcell 1, adapter with tag ATCACG.
|   |-- SN_11_LPS: Sample from person 2, treated with LPS
|   |   `-- Paired end reads in fastq format from four lanes of flowcell 1, adapter with tag TGACCA.
|   |-- SN_11_UNST: Sample from person 2, untreated
|   |   `-- Paired end reads, fastq format from four lanes of flowcell 1, adapter with tag TTAGGC.
|   |-- SN_12_LPS: Sample from person 3, treated with LPS
|   |   `-- Paired end reads in fastq format from four lanes of flowcell 1, adapter with tag GCCAAT.
|   |-- SN_12_UNST: Sample from person 3, untreated
|   |   `-- Paired end reads, fastq format from four lanes of flowcell 1, adapter with tag ACAGTG.
`-- Flowcell 2: 130104_SN866_0198_BC1DAYACXX
    |-- SN_10_LPS: Sample from person 1, treated with LPS
    |   `-- Paired end reads in fastq format from four lanes of flowcell 2, adapter with tag CGATGT.
    |-- SN_10_UNST: Sample from person 1, untreated
    |   `-- Paired end reads in fastq format from four lanes of flowcell 2, adapter with tag ATCACG.
    |-- SN_11_LPS: Sample from person 2, treated with LPS
    |   `-- Paired end reads in fastq format from four lanes of flowcell 2, adapter with tag TGACCA.
    |-- SN_11_UNST: Sample from person 2, untreated
    |   `-- Paired end reads, fastq format from four lanes of flowcell 2, adapter with tag TTAGGC.
    |-- SN_12_LPS: Sample from person 3, treated with LPS
    |   `-- Paired end reads in fastq format from four lanes of flowcell 2, adapter with tag GCCAAT.
    |-- SN_12_UNST: Sample from person 3, untreated
        `-- Paired end reads, fastq format from four lanes of flowcell 2, adapter with tag ACAGTG.
```

####Project plan
We came up with a plan for the first part of the experiment, described below.

#####Subsampling

After seeing how the data is structured, we decided that the easiest way to do subsampling would be to select only a subset of the lanes. By first using lane 1 from both flowcells, then using lane 1 and 2, then 1,2 and 3 and so forth. This would let us compare the results of using 25%, 50%, 75% and 100% of the data.

#####Preprocessing

First, we'll trim possible adaptor sequences from the reads. Often, the adaptors are sequenced partially or completely. Keeping them of the reads complicates mapping, and hence they should be removed. [Cutadapt](https://cutadapt.readthedocs.org/en/stable/) is a commonly used tool to remove adapter sequences, primers, poly-A tails and other types of unwanted artifacts from reads. Cutadapt can take the adapter tag sequence stated in the filename of the fastq files as input. Low quality reads should also be removed, and [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is commonly used filter reads based on quality. FastQC also provides analysis of reads based on quality and other factors. A [preliminary analysis of `SN10_UNST_ATCACG_L001_R2`](https://cdn.rawgit.com/aerugo/BB2490-RNASeq-Project/master/results/2016-02-23/20160223-SN10_UNST_ATCACG_L001_R2-fastqresult.html) with FastQC shows that avarage base quality score per read is over 25 for the vast majority of reads, and that there is a non-random distribution of per base sequence content at the beginning of many of the reads. This indicates that there is probably an adapter sequence present. 

The wrapper [TrimGalore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) combines cutadapt and fastQC in a single command. We used this in a [bashscript](https://github.com/aerugo/BB2490-RNASeq-Project/blob/master/src/trimming.sh), which we tried out on `SN10_UNST_ATCACG_L001_R2` to compare with the FastQC analysis of the raw reads. We put this job in queue on Uppmax.

#####Mapping and post-processing

We are planning to use STAR for mapping and will then filter based on mapping quality with SAMTools.

#####Differential Expression

After the steps above are complete, we will have 6 BAM files - treated and untreated samples from three subjects.

#####Allele specific expression

How this will be done is yet to be determined.
