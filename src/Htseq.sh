#! /bin/bash -l
#SBATCH -A g2015056
#SBATCH -t 40:00:00
#SBATCH -J Count_reads
#SBATCH -p core -n 4
#SBATCH -C thin
#SBATCH -e /proj/g2015056/nobackup/Subsample/Htseq/HTseq_errSN10.txt
#SBATCH -o /proj/g2015056/nobackup/Subsample/Htseq/HTSeq_outSN10.txt
#SBATCH --mail-type=All
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se


## Count reads in each features

module load bioinfo-tools
module add htseq

input=/proj/g2015056/BB2490/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/20160304-starmapping/
output=/proj/g2015056/nobackup/Subsample/Htseq/

SAMPLES=(SN10_UNST_ATCACG_ SN_10_LPS_CGATGT SN11_UNST_TTAGGC SN_11_LPS_TGACCA SN_12_LPS_GCCAAT SN12_UNST_ACAGTG)

cd ${input}

for i in "${SAMPLES[@]}"; do
    for file in ${i}*"out.sam"; do
        echo $file
        substring=${file:0:17}
        name=${substring}"count.txt"
        echo $name
        htseq-count -m intersection-strict -t exon -i gene_id ${file} /proj/g2015056/nobackup/Homo_sapiens.GRCh38.83.gtf > ${output}${name} &
    done
    wait
done