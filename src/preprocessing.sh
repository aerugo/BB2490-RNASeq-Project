#!/bin/bash

# This script runs preprocessing on fastq data and saves it in results

module load bioinfo-tools
module add FastQC cutadapt python TrimGalore

trim_galore -a CGATGT -q 15 -s 5 -e 0.05 --length 48 <fastq_file>