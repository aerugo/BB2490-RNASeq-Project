#! /bin/bash -l
#SBATCH -A g2015056
#SBATCH -t 4:00:00
#SBATCH -J trim and combine
#SBATCH -p node -n 8
#SBATCH -C thin
#SBATCH -e /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/trim_all_err.txt
#SBATCH -o /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/trim_all_out.txt
#SBATCH --mail-type=All
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se


## Trimming

input=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/20160302/combined/
output=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/20160302/trimmed/

cd ${input}
for file in *.fastq.gz; do
    echo "first file"
    filename=${basename ${file}}
    substring=${filename:0:22}
    echo ${substring}"R1_001_val_1.fastq"
    ls ${substring}"R1_001_val_1.fastq"
    #trim_galore -q 25 --stringency 1 --paired --length 25 -o ${output} ${read1} ${read2}
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