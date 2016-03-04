#! /usr/bin/env python
# This script takes STAR logfiles and retrieves the number of mapped
# reads. Some not-so-elegant solutions are marked with "### UNSTABLE 
# SOLUTION"
# Syntax: ./counter.py [DIRECTORY_FOR_INPUT_BAM/LOG] [OUTPUT_FILENAME]

import os
import sys

# Retrieves a file list containing all the log and bam files in
# the chosen directory. After that, the list is sorted.
file_list = []
directory = os.listdir(sys.argv[1])
searchstring1 = 'Log.final.out'
searchstring2 = 'Aligned.out.bam'
for fname in directory:
    if searchstring1 in fname or searchstring2 in fname:
        file_list.append(fname)
file_list.sort()

# Devides file list into bam and log files.
file_list_bam = file_list[
                0:len(file_list):2]  ### UNSTABLE SOLUTION, depends the list ordering to be [bam, log, bam, log...]
file_list_log = file_list[
                1:len(file_list):2]  ### UNSTABLE SOLUTION, depends the list ordering to be [bam, log, bam, log...]

# Opens an output file.
with open(sys.argv[2], 'w') as out:
    previous_name = ''
    counter = -1  ### UNSTABLE SOLUTION, depends on list order

    # Loop for every log file found.
    for file_name in file_list_bam:
        counter = counter + 1  ### UNSTABLE SOLUTION, see line 29 (counter)

        # Writes a header into the output containing the sample
        # name and the treatment. If the sample and treatment are
        # the same as they were previous loop it skips this step.
        if file_name.split('_')[:2] == previous_name:
            pass
        else:
            out.writelines('\n#')
            for elements in file_list_bam[counter].split('_')[:2]:  ### UNSTABLE SOLUTION, see line 29 (counter)
                out.writelines(str(elements) + ' ')

        # Opens the file and retrieves the number of mapped reads
        # and writes to output file.
        count_column = open(file_name, 'r').readlines()[
            8]  ### UNSTABLE SOLUTION, takes 8th column and does not look for "Uniquely mapped read count"
        out.writelines('\n' + file_list_bam[counter] + '\t' + count_column.split()[
            -1])  ### UNSTABLE SOLUTION, see line 29 (counter)

        # Updates previous name variable.
        previous_name = file_name.split('_')[:2]
