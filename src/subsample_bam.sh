#! /bin/bash -l
#SBATCH -A g2015056
#SBATCH -t 40:00:00
#SBATCH -J subsample
#SBATCH -p core -n 4
#SBATCH -C thin
#SBATCH -e /proj/g2015056/nobackup/Subsample/subsample_err.txt
#SBATCH -o /proj/g2015056/nobackup/Subsample/subsample_out.txt
#SBATCH --mail-type=All
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se

#module load bioinfo-tools
#module add samtools

#output=/proj/g2015056/nobackup/Subsample/
input=home/hugia/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/20160304-starmapping
output_subsample=${output}
output=/Users/aerugo

subsample_status=$(python generate_subsampling.py)

if [ ${subsample_status} = "success" ]
then
    echo "Do you want to generate subsample files? Choose list to see overview."
    select yn in "Yes" "No" "List"; do
        case $yn in
            Yes ) :; break;;
            No ) exit;;
            List ) less -FX subsample_scheme.txt; echo; echo "Do you want to generate subsample files? Choose list to see overview." ;;
        esac
    done
fi

now=$(date +"%m_%d_%Y")
mkdir -p ${output}/${now}

while read p; do
  if [[ $p == @subsample* ]]
  then
    mkdir -p ${output}/${now}/${p#?}
    output_subsample=${output}/${now}/${p#?}
  else
    IFS=' ' read -ra line <<< "$p"
    name=${line[0]}
    original_count=${line[1]}
    sample_count=${line[2]}
    if [[ ${original_count} == ${sample_count} ]]
    then
        fraction=$(bc <<< "scale=4;$sample_count/$original_count")
        echo "ln -s ${input}/infile.bam > ${output_subsample}/outfile.bam"
    else
        fraction=$(bc <<< "scale=4;$sample_count/$original_count")
        echo "samtools view -b -s ${fraction} ${input}/infile.bam > ${output_subsample}/outfile.bam"
    fi
  fi
done <subsample_scheme.txt