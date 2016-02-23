#! /bin/bash -l
#SBATCH -A g2015056
#SBATCH -t 2:00:00
#SBATCH -J Trim_test
#SBATCH -p node -n 1
#SBATCH -C thin
#SBATCH -e /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/trim_err.txt
#SBATCH -o /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/trim_out.txt
#SBATCH --mail-type=All
#SBATCH --mail-user=hugia@kth.se

module add bioinfo-tools
module add python/2.7
module load FastQC

output=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/
input=/proj/g2015056/BB2490/proj5_ASE/130104_SN866_0197_AC1DLVACXX/Sample_SN10_UNST/

#####To check: Better to use those?
## nextera1=CTGTCTCTTATACACATCTGACGCTGCCGACGA
## nextera2=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
#####

adapter=ATCACG

##Nextera
mkdir /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/fastqc_results
fastqc -o /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/fastqc_results -f /proj/g2015056/BB2490/proj5_ASE/130104_SN866_0197_AC1DLVACXX/Sample_SN10_UNST/SN10_UNST_ATCACG_L001_R1_001.fastq.gz
fastqc -o /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/fastqc_results -f /proj/g2015056/BB2490/proj5_ASE/130104_SN866_0197_AC1DLVACXX/Sample_SN10_UNST/SN10_UNST_ATCACG_L001_R1_002.fastq.gz
wait

