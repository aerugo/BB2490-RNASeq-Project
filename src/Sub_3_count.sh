#! /bin/bash -l
#SBATCH -A g2015056
#SBATCH -t 24:00:00
#SBATCH -J HTseq_subsample_3
#SBATCH -p core -n 6
#SBATCH -C thin
#SBATCH -e /proj/g2015056/nobackup/Subsample/Htseq/HTseq_Subsample_3_err.txt
#SBATCH -o /proj/g2015056/nobackup/Subsample/Htseq/HTSeq_Subsample_3_out.txt
#SBATCH --mail-type=All
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se

## HTseqcount
## Count reads in each features

module load bioinfo-tools
module add htseq

input=/proj/g2015056/nobackup/Subsample/03_10_2016/subsample_3/
output=/proj/g2015056/nobackup/Subsample/Htseq/Subsample_3/

SAMPLES=(SN_11_LPS_TGACCA)

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
