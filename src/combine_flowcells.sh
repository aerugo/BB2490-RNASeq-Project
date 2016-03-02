#! /bin/bash -l                                                                                                                                                                                                                              
#SBATCH -A g2015056                                                                                                                                                                                     
#SBATCH -t 4:00:00                                                                                                                                                                                                                           
#SBATCH -J trim and combine
#SBATCH -p node -n 8
#SBATCH -C thin                                                                                                                                                                                                                              
#SBATCH -e /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/combine_flowcells_err.txt
#SBATCH -o /proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/combine_flowcells_out.txt
#SBATCH --mail-type=All
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se

# Remember to make directory

module load bioinfo-tools
module add python/2.7
module add FastQC
module add cutadapt
module add TrimGalore

# First, combine all same lanes from flow cells into into one data set each

output=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/preprocessed_reads/20160302/combined/
input_fc1=/proj/g2015056/BB2490/proj5_ASE/130104_SN866_0197_AC1DLVACXX/
input_fc2=/proj/g2015056/BB2490/proj5_ASE/130104_SN866_0198_BC1DAYACXX/

## Combine with cat
# Sample_SN_10_LPS

cd ${output}

for directory in ${input_fc1}/S*; do
    cd ${directory}
    for file in *.fastq.gz; do
        filename=$(basename "$file")
        echo "File 1"
        echo ${filename}
        file2=${input_fc2}$(basename "$directory")$(basename "$file")
        echo ${file2}
        ls ${file2}
        #echo "Yay1"
        #ls ${input_fc1}/${directory##*/}/${filename}
        #echo "Yay2"
        #echo ${input_fc2}/${directory##*/}/${filename}
        #echo "Yay3"
        #ls ${input_fc2}/${directory##*/}/${filename}
        #echo "Yay4"
    done
done

# cat ${input_fc1}/Sample_SN_10_LPS/SN_10_LPS_CGATGT_L001_R1_001.fastq.gz ${input_fc2}/Sample_SN_10_LPS/SN_10_LPS_CGATGT_L001_R1_001.fastq.gz
