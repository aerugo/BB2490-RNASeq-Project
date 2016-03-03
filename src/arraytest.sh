#! /bin/bash -l
#SBATCH -A g2015056
#SBATCH -t 4:00:00
#SBATCH -J trimall
#SBATCH -p core -n 2
#SBATCH -C thin
#SBATCH -e /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/trim_all_err.txt
#SBATCH -o /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/trim_all_out.txt
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

SAMPLES=(SN_11_LPS_TGACCA SN11_UNST_TTAGGC SN_12_LPS_GCCAAT SN12_UNST_ACAGTG)

cd ${input}

for i in "${SAMPLES[@]}"; do
    for file in ${i}*"R1_001.fastq"; do
        substring=${file:0:22}
        paired_read1=${substring}"R1_001.fastq"
        paired_read2=${substring}"R2_001.fastq"
        echo "trim_galore -q 25 --stringency 1 --paired --length 25 --fastqc -o "${output}" "${paired_read1}" "${paired_read2} &
    done
    wait
done

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
