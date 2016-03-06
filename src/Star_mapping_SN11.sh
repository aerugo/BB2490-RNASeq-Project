#! /bin/bash -l                                                                                                                                                                                                                              
#SBATCH -A b2012058
#SBATCH -t 8:00:00                                                                                                                                                                                                                           
#SBATCH -J bamtosam
#SBATCH -p node -n 8
#SBATCH -C thin                                                                                                                                                                                                                              
#SBATCH -e /proj/b2012058/nobackup/private/Exome_data/Whole_genome_sequuencing/scripts/info/samtobam_err.txt
#SBATCH -o /proj/b2012058/nobackup/private/Exome_data/Whole_genome_sequuencing/scripts/infoout.txt
#SBATCH --mail-type=All
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se

module load bioinfo-tools star

#mkdir -p /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/20160304-starmapping
reference=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_index/20160226-starindex/
input=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/20160302/trimmed/
output=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/20160304-starmapping/

##Script
STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN_11_LPS_TGACCA_L002_R1_001_val_1.fq ${input}SN_11_LPS_TGACCA_L002_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN_11_LPS_TGACCA_L002

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN_11_LPS_TGACCA_L003_R1_001_val_1.fq ${input}SN_11_LPS_TGACCA_L003_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN_11_LPS_TGACCA_L003

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN_11_LPS_TGACCA_L004_R1_001_val_1.fq ${input}SN_11_LPS_TGACCA_L004_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN_11_LPS_TGACCA_L004

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN_11_LPS_TGACCA_L005_R1_001_val_1.fq ${input}SN_11_LPS_TGACCA_L005_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN_11_LPS_TGACCA_L005


STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN11_UNST_TTAGGC_L002_R1_001_val_1.fq ${input}SN11_UNST_TTAGGC_L002_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN11_UNST_TTAGGC_L002

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN11_UNST_TTAGGC_L003_R1_001_val_1.fq ${input}SN11_UNST_TTAGGC_L003_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN11_UNST_TTAGGC_L003

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN11_UNST_TTAGGC_L003_R1_001_val_1.fq ${input}SN11_UNST_TTAGGC_L004_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN11_UNST_TTAGGC_L004

STAR --genomeDir ${reference} --runThreadN 8 --readFilesIn ${input}SN11_UNST_TTAGGC_L003_R1_001_val_1.fq ${input}SN11_UNST_TTAGGC_L005_R2_001_val_2.fq --outReadsUnmapped Fastx --outFileNamePrefix ${output}/SN11_UNST_TTAGGC_L005