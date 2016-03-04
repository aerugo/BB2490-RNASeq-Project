#! /bin/bash -l                                                                                                                                                                                                                              
#SBATCH -A g2015056                                                                                                                                                                                     
#SBATCH -t 8:00:00                                                                                                                                                                                                                           
#SBATCH -J star_mapping-S12
#SBATCH -p node -n 8
#SBATCH -C thin                                                                                                                                                                                                                              
#SBATCH -e /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/star_mapping_Mar_02err.txt
#SBATCH -o /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/star_mapping_Mar_02_out.txt
#SBATCH --mail-type=All
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se

module load bioinfo-tools star

#mkdir -p /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/20160304-starmapping
reference=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_index/20160226-starindex/
input=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/20160302/trimmed/
output=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/20160304-starmapping/

##Script
STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN_12_LPS_GCCAAT_L003_R1_001_val_1.fq ${input}SN_12_LPS_GCCAAT_L003_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN_12_LPS_GCCAAT_L003

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN_12_LPS_GCCAAT_L004_R1_001_val_1.fq ${input}SN_12_LPS_GCCAAT_L004_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN_12_LPS_GCCAAT_L004

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN_12_LPS_GCCAAT_L005_R1_001_val_1.fq ${input}SN_12_LPS_GCCAAT_L005_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN_12_LPS_GCCAAT_L005

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN_12_LPS_GCCAAT_L006_R1_001_val_1.fq ${input}SN_12_LPS_GCCAAT_L006_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN_12_LPS_GCCAAT_L006


STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN12_UNST_ACAGTG_L003_R1_001_val_1.fq ${input}SN12_UNST_ACAGTG_L003_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN12_UNST_ACAGTG_L003

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN12_UNST_ACAGTG_L004_R1_001_val_1.fq ${input}SN12_UNST_ACAGTG_L004_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN12_UNST_ACAGTG_L004
	
STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN12_UNST_ACAGTG_L005_R1_001_val_1.fq ${input}SN12_UNST_ACAGTG_L005_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN12_UNST_ACAGTG_L005

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN12_UNST_ACAGTG_L006_R1_001_val_1.fq ${input}SN12_UNST_ACAGTG_L006_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN12_UNST_ACAGTG_L006
