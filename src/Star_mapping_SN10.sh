#! /bin/bash -l                                                                                                                                                                                                                              
#SBATCH -A g2015056                                                                                                                                                                                     
#SBATCH -t 6:00:00
#SBATCH -J star_mapping-S10
#SBATCH -p node -n 8
#SBATCH -C thin                                                                                                                                                                                                                              
#SBATCH -e /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/star_mapping_Mar_04_err.txt
#SBATCH -o /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/star_mapping_Mar_04_out.txt
#SBATCH --mail-type=All
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se

module load bioinfo-tools star

mkdir -p /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/20160304-starmapping
reference=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_index/20160226-starindex/
input=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/20160302/trimmed/
output=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/20160304-starmapping/

##Script

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN10_UNST_ATCACG_L001_R1_001_val_1.fq ${input}SN10_UNST_ATCACG_L001_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN10_UNST_ATCACG_L001

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN10_UNST_ATCACG_L002_R1_001_val_1.fq ${input}SN10_UNST_ATCACG_L002_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN10_UNST_ATCACG_L002

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN10_UNST_ATCACG_L003_R1_001_val_1.fq ${input}SN10_UNST_ATCACG_L003_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN10_UNST_ATCACG_L003

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN10_UNST_ATCACG_L004_R1_001_val_1.fq ${input}SN10_UNST_ATCACG_L004_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN10_UNST_ATCACG_L004


STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN_10_LPS_CGATGT_L001_R1_001_val_1.fq ${input}SN_10_LPS_CGATGT_L001_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN_10_LPS_CGATGT_L001

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN_10_LPS_CGATGT_L001_R1_001_val_1.fq ${input}SN_10_LPS_CGATGT_L002_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN_10_LPS_CGATGT_L002

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN_10_LPS_CGATGT_L001_R1_001_val_1.fq ${input}SN_10_LPS_CGATGT_L003_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN_10_LPS_CGATGT_L003

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN_10_LPS_CGATGT_L001_R1_001_val_1.fq ${input}SN_10_LPS_CGATGT_L004_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN_10_LPS_CGATGT_L004 
