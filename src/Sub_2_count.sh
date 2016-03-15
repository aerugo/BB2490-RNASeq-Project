#! /bin/bash -l
#SBATCH -A g2015056
#SBATCH -t 24:00:00
#SBATCH -J HTseq_subsample_2
#SBATCH -p core -n 6
#SBATCH -C thin
#SBATCH -e /proj/g2015056/nobackup/Subsample/Htseq/HTseq_Subsample_2_err.txt
#SBATCH -o /proj/g2015056/nobackup/Subsample/Htseq/HTSeq_Subsample_2_out.txt
#SBATCH --mail-type=All
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se

## HTseqcount
## Count reads in each features

module load bioinfo-tools
module add htseq

input=/proj/g2015056/nobackup/Subsample/03_10_2016/subsample_2/
output=/proj/g2015056/nobackup/Subsample/Htseq/Subsample_2/

SAMPLES=(SN11_UNST_TTAGGC SN_11_LPS_TGACCA SN_12_LPS_GCCAAT SN12_UNST_ACAGTG)

cd ${input}

for i in "${SAMPLES[@]}"; do
    for file in ${i}*"out.sam"; do
        echo $file
        substring=${file:0:21}
        name=${substring}"_count.txt"
        echo $name
	htseq-count -f bam -m intersection-strict -t exon -i gene_id ${file}  /proj/g2015056/nobackup/Homo_sapiens.GRCh38.83.gtf  > ${output}${name} &
    done
    wait
done
