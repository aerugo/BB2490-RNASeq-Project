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

# ------------------------
# This script generates subsamples of a set of bam files in the input directory.
# Subsamples and a log file are generated in the output directory.
# Files to be processed must be listed in a separate file formatted as follows
#               # Subsample_name_1
#               # filename_1.bam    number_of_reads
#               # filename_2.bam    number_of_reads
#               ...
#               # Subsample_name_2
#               # filename_1.bam    number_of_reads
#               # filename_2.bam    number_of_reads
#               ...
# Settings and path to list of bam file given in config.ini, which must be present in same file as python script.
# ------------------------

module load bioinfo-tools
module add samtools

output=/proj/g2015056/nobackup/Subsample/
input=/proj/g2015056/proj5_ASE/BB2490-RNASeq-Project/data/star_mapping/20160304-starmapping
output_subsample=${output}

rm subsample_scheme.txt
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
else
    exit
fi

now=$(date +"%m_%d_%Y")
mkdir -p ${output}/${now}

touch tempfile.txt
printf "\n \n" >> tempfile.txt
printf "Subsampling log\n" >> tempfile.txt

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
        echo "ln -s ${input}/${name} > ${output_subsample}/${name}"
        printf "ln -s ${input}/${name} > ${output_subsample}/${name} \n" >> tempfile.txt # Should print command to tempfile
    else
        fraction=$(bc <<< "scale=4;$sample_count/$original_count")
        echo "samtools view -b -s 0${fraction} ${input}/${name} > ${output_subsample}/${name}"
        printf "samtools view -b -s 0${fraction} ${input}/${name} > ${output_subsample}/${name} \n" >> tempfile.txt
    fi
  fi
done <subsample_scheme.txt

cat subsample_scheme.txt tempfile.txt > ${output}/${now}/subsample_log.txt
rm tempfile.txt