#! /bin/bash -l
#SBATCH -A g2015056
#SBATCH -t 4:00:00
#SBATCH -J star_indexing
#SBATCH -p node -n 8 
#SBATCH -C thin 
#SBATCH -e /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_index/star_indexing_err.txt
#SBATCH -o /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_index/star_indexing_out.txt
#SBATCH --mail-type=All
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se 

module load bioinfo-tools star

mkdir -p /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_index/20160226-starindex
output=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_index/20160226-starindex
input=/proj/g2015056/nobackup/

##Script
STAR --runMode genomeGenerate --genomeDir ${output} --genomeFastaFiles ${input}Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile ${input}Homo_sapiens.GRCh38.83.gtf --runThreadN 8
