#! /bin/bash -l
#SBATCH -A g2015056
#SBATCH -t 10:00:00
#SBATCH -J trimingandfastq_sample10
#SBATCH -p core -n 8
#SBATCH -C thin
#SBATCH -e /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/20160302/trimmed/trim_all_err.txt
#SBATCH -o /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/20160302/trimmed/trim_out_err.txt
#SBATCH --mail-type=All
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se

## Trimming

module load bioinfo-tools
module add python/2.7
module add FastQC
module add cutadapt
module add TrimGalore

input=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/20160302/combined/
output=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/20160302/trimmed/


trim_galore -q 25 --stringency 1 --paired --length 25 --fastqc -o ${output} ${input}SN11_UNST_TTAGGC_L002_R1_001.fastq ${input}SN11_UNST_TTAGGC_L002_R2_001.fastq &

trim_galore -q 25 --stringency 1 --paired --length 25 --fastqc -o ${output} ${input}SN11_UNST_TTAGGC_L003_R1_001.fastq ${input}SN11_UNST_TTAGGC_L003_R2_001.fastq &

trim_galore -q 25 --stringency 1 --paired --length 25 --fastqc -o ${output} ${input}SN11_UNST_TTAGGC_L004_R1_001.fastq ${input}SN11_UNST_TTAGGC_L004_R2_001.fastq &

trim_galore -q 25 --stringency 1 --paired --length 25 --fastqc -o ${output} ${input}SN11_UNST_TTAGGC_L005_R1_001.fastq ${input}SN11_UNST_TTAGGC_L005_R2_001.fastq &
wait


trim_galore -q 25 --stringency 1 --paired --length 25 --fastqc -o ${output} ${input}SN_11_LPS_TGACCA_L002_R1_001.fastq ${input}SN_11_LPS_TGACCA_L002_R1_001.fastqt &

trim_galore -q 25 --stringency 1 --paired --length 25 --fastqc -o ${output} ${input}SN11_UNST_ATCACG_L002_R1_001.fastq ${input}SN11_UNST_ATCACG_L002_R2_001.fastq &

trim_galore -q 25 --stringency 1 --paired --length 25 --fastqc -o ${output} ${input}SN11_UNST_ATCACG_L003_R1_001.fastq ${input}SN11_UNST_ATCACG_L003_R2_001.fastq &

trim_galore -q 25 --stringency 1 --paired --length 25 --fastqc -o ${output} ${input}SN11_UNST_ATCACG_L004_R1_001.fastq ${input}SN11_UNST_ATCACG_L004_R2_001.fastq &
wait



#module load bioinfo-tools FastQC
#fastqc ${output}SN_10_LPS_CGATGT_L001_R1_001_val_1.fq.gz &
#fastqc ${output}SN_10_LPS_CGATGT_L001_R2_001_val_2.fq.gz &
#wait

#fastqc ${output}SN10_UNST_ATCACG_L001_R1_001_val_1.fq.gz &
#fastqc ${output}SN10_UNST_ATCACG_L001_R2_001_val_2.fq.gz &
#wait
#
#mkdir -p /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/results/2016-02-26/fastq_results/
#mv -v ${output}*.zip /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/results/2016-02-26/fastq_results/
#mv -v ${output}*.html /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/results/2016-02-26/fastq_results/