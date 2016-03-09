#! /usr/bin/env python3
import os
import sys


# Parse arguments and get config data


def read_config():
    s = open(os.path.join(__location__, './config.ini'), 'r').read()
    return eval(s)


# Clustering is calculated by first generating a dictionary of the bam files and then applying the algorithm
# Don't worry about bam_dict - file_list is what we are using!

def generate_bam_dict():
    bam_dict = {}
    file_list = []
    with open(bam_list) as f:
        sample = ""
        read_counts = []
        for line in f:
            line = line.strip()
            if str.strip(line) != "":
                if line[0] == "#":
                    sample = line[1:]
                    bam_dict[sample] = {}
                    bam_dict[sample]["file_list"] = []
                    read_counts = []
                else:
                    file_name = str.split(line)[0]
                    read_count = int(float(str.split(line)[1]))
                    read_counts.append(read_count)
                    file_list.append([file_name, read_count, sample])
                    bam_dict[sample]["file_list"].append([file_name, read_count])
                    bam_dict[sample]["average_read_count"] = sum(read_counts) / len(read_counts)
                    bam_dict[sample]["lowest_read_count"] = min(read_counts)
                    bam_dict[sample]["highest_read_count"] = max(read_counts)

    return bam_dict, file_list


def calculate_clustering(bam_dict, cut):
    cut *= 1000000
    tuples = []
    cut_list = []
    max_read_count = 0
    for i, k in bam_dict.items():
        tuples.append((k["average_read_count"], k["highest_read_count"]))
    tuples.sort(key=lambda tup: tup[0])
    last_average = tuples[0][0]
    for i in tuples:
        if i[0] - last_average > cut:
            cut_list.append(max_read_count)
        last_average = i[0]
        max_read_count = i[1]
    cut_list.append(max_read_count)
    return cut_list


# Outputs an array of arrays. Each array is a subsample scheme.
# Each subsample scheme array contains one tuple per bam file.
# Each tuple contains (filename, original read count, sub sample read count)


def return_subsample_data(smallest_cut):
    bam_dict, file_list = generate_bam_dict()
    cut_list = calculate_clustering(bam_dict, smallest_cut)
    subsample_data = []

    output = open(os.path.join(__location__, './subsample_scheme.txt'), 'w')

    for i, cut in enumerate(cut_list):
        output.write("@" + "subsample_" + str(i+1) + "\n")
        subsample = []
        for file in file_list:
            if file[1] <= cut:
                bam_file = [file[0], file[1], file[1]]
                output.write(bam_file[0] + " " + str(bam_file[1]) + " " + str(bam_file[2]) + " " + str(file[2]) + "\n")
            else:
                bam_file = [file[0], file[1], cut]
                output.write(bam_file[0] + " " + str(bam_file[1]) + " " + str(bam_file[2]) + " " + str(file[2]) + "\n")
            subsample.append(bam_file)
        subsample_data.append(subsample)

    return True


if __name__ == "__main__":

    __location__ = os.path.realpath(
        os.path.join(os.getcwd(), os.path.dirname(__file__)))

    # Loading config
    CONFIG = read_config()
    cut = int(CONFIG['cut'])
    bam_list = str(os.path.join(__location__, CONFIG['bam_list']))
    smallest_fraction = CONFIG['smallest_fraction']

    data = return_subsample_data(cut)
    if data:
        sys.stdout.write("success")
