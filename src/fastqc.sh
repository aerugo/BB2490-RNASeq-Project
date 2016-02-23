#! /bin/bash -l

module add bioinfo-tools
module add python
module add FastQC
module add cutadapt
module add TrimGalore

output=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/
input=/proj/g2015056/BB2490/proj5_ASE/130104_SN866_0197_AC1DLVACXX/Sample_SN10_UNST/

#####To check: Better to use those?
## nextera1=CTGTCTCTTATACACATCTGACGCTGCCGACGA
## nextera2=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
#####

# mkdir /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/results/fastqc_results
fastqc -o /home/hugia/proj5_ASE/BB2490-RNASeq-Project/results/2016-02-23/fastq_results_trimgalore_SN10_UNST_ATCACG_L001_R2 /home/hugia/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/SN10_UNST_ATCACG_L001_R1_001_val_1.fq.gz
#fastqc -o /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/results/fastqc_results /proj/g2015056/BB2490/proj5_ASE/130104_SN866_0197_AC1DLVACXX/Sample_SN10_UNST/SN10_UNST_ATCACG_L001_R2_001.fastq.gz
