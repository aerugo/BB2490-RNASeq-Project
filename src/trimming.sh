#! /bin/bash -l                                                                                                                                                                                                                              
#SBATCH -A g2015056                                                                                                                                                                                     
#SBATCH -t 4:00:00                                                                                                                                                                                                                           
#SBATCH -J Preprocessing                                                                                                                                                                                                                    
#SBATCH -p core -n 2                                                                                                                                                                                                                        
#SBATCH -C thin                                                                                                                                                                                                                              
#SBATCH -e /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/trim_err_2.txt                                                                                                                             
#SBATCH -o /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/trim_out_2.txt                                
#SBATCH --mail-type=All                                                                                                                                                                                                                      #SBATCH --mail-user=sailendra.pradhananga@scilifelab.se                                                                                                                                                                                      

module load bioinfo-tools
module add python/2.7
module add FastQC
module add cutadapt
module add TrimGalore

mkdir -p /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/20160226-trimgalore
output=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/20160226-trimgalore/
input=/proj/g2015056/BB2490/proj5_ASE/130104_SN866_0197_AC1DLVACXX/

#####To check: Better to use those?                                                                                                                                                                                                           
## nextera1=CTGTCTCTTATACACATCTGACGCTGCCGACGA                                                                                                                                                                                                 
## nextera2=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC                                                                                                                                                                                                
#####                                                                                                                                                                                                                                         

adapter=ATCACG

##Nextera                                                                                                                                                                                                                                     
trim_galore -q 25 --stringency 1 --paired --length 25 -o ${output} ${input}/Sample_SN_10_LPS/SN_10_LPS_CGATGT_L001_R1_001.fastq.gz ${input}/Sample_SN_10_LPS/SN_10_LPS_CGATGT_L001_R2_001.fastq.gz &
trim_galore -q 25 --stringency 1 --paired --length 25 -o ${output} ${input}/Sample_SN10_UNST/SN10_UNST_ATCACG_L001_R1_001.fastq.gz ${input}/Sample_SN10_UNST/SN10_UNST_ATCACG_L001_R2_001.fastq.gz &
wait

module load bioinfo-tools FastQC
fastqc ${output}SN_10_LPS_CGATGT_L001_R1_001_val_1.fq.gz &
fastqc ${output}SN_10_LPS_CGATGT_L001_R2_001_val_2.fq.gz &
wait

fastqc ${output}SN10_UNST_ATCACG_L001_R1_001_val_1.fq.gz &
fastqc ${output}SN10_UNST_ATCACG_L001_R2_001_val_2.fq.gz &
wait

mkdir -p /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/results/2016-02-26/fastq_results/
mv -v ${output}*.fq_fastqc* /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/results/2016-02-26/fastq_results/

#trim_galore -a ${adapter} -a2 ${adapter} -q 25 --stringency 5 --paired --length 25 -o ${output} ${input}S0156_1.fastq.gz ${input}S0156_2.fastq.gz &
#trim_galore -a ${adapter} -a2 ${adapter} -q 25 --stringency 5 --paired --length 25 -o ${output} ${input}S0160_1.fastq.gz ${input}S0160_2.fastq.gz &
#wait

