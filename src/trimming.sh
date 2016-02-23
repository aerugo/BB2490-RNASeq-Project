#! /bin/bash -l                                                                                                                                                                                                                               
#SBATCH -A g2015056                                                                                                                                                                                                                           
#SBATCH -t 2:00:00                                                                                                                                                                                                                           
#SBATCH -J Trim_test                                                                                                                                                                                                                         
#SBATCH -p node -n 1                                                                                                                                                                                                                          
#SBATCH -C thin                                                                                                                                                                                                                               
#SBATCH -e /proj/g2015056/nobackup/info/trim_err.txt                                                                                                                                                                      
#SBATCH -o /proj/b2012056/nobackup/info/trim_out.txt                                                                                                                                                                      
#SBATCH --mail-type=All                                                                                                                                                                                                                       
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se                                                                                                                                                                                      

module add bioinfo-tools
module add python/2.7
module load FastQC
module load cutadapt
module load TrimGalore


#####To check: Better to use those?                                                                                                                                                                                                           
## nextera1=CTGTCTCTTATACACATCTGACGCTGCCGACGA                                                                                                                                                                                                 
## nextera2=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC                                                                                                                                                                                                
#####                                                                                                                                                                                                                                         

nextera=CTGTCTCTTATACACATCT

##Nextera                                                                                                                                                                                                                                     
trim_galore -a ${nextera} -a2 ${nextera} -q 25 --stringency 5 --paired --length 25 -o ${output} ${input}S0143_1.fastq.gz ${input}S0143_2.fastq.gz &
trim_galore -a ${nextera} -a2 ${nextera} -q 25 --stringency 5 --paired --length 25 -o ${output} ${input}S0156_1.fastq.gz ${input}S0156_2.fastq.gz &
trim_galore -a ${nextera} -a2 ${nextera} -q 25 --stringency 5 --paired --length 25 -o ${output} ${input}S0160_1.fastq.gz ${input}S0160_2.fastq.gz &
wait

