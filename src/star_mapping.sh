#! /bin/bash -l                                                                                                                                                                                                                              
#SBATCH -A g2015056                                                                                                                                                                                     
#SBATCH -t 4:00:00                                                                                                                                                                                                                           
#SBATCH -J star_indexing
#SBATCH -p node -n 8
#SBATCH -C thin                                                                                                                                                                                                                              
#SBATCH -e /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/star_indexing_err.txt
# SBATCH -o /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/star_indexing_out.txt
# SBATCH --mail-type=All
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se

module load bioinfo-tools
module add python/2.7
module add star

mkdir -p /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/20160226-starmapping
reference=/home/hugia/proj5_ASE/BB2490-RNASeq-Project/data/reference_genome_indexes/20160226-starindex/
input=/home/hugia/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/20160226-trimgalore/
output=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/20160226-starmapping/

##Script
STAR --genomeDir ${reference}CHANGETOINDEX --runThreadN 8 --readFilesIn ${input}SN10_UNST_ATCACG_L001_R1_001_val_1.fq.gz ${input}SN10_UNST_ATCACG_L001_R2_001_val_2.fq.gz --readFilesCommand gunzip -c --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN10_UNST_ATCACG_L001/
STAR --genomeDir ${reference}CHANGETOINDEX --runThreadN 8 --readFilesIn ${input}SN_10_LPS_CGATGT_L001_R1_001_val_1.fq.gz ${input}SN_10_LPS_CGATGT_L001_R2_001_val_2.fq.gz --readFilesCommand gunzip -c --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN_10_LPS_CGATGT_L001/