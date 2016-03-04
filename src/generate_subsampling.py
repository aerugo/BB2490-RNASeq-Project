#! /usr/bin/env python3

import sys
from operator import itemgetter

# Parse arguments and get config data


def read_config():
    s = open('config.ini', 'r').read()
    return eval(s)


def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cut", help="Difference in average number of reads between samples required to cluster.")
    parser.add_argument("-d", "--directory", help="Directory containing bam to subsample")
    parser.add_argument("-l", "--bam_list", help="Text file containing read counts for each bam file to be divided")
    parser.add_argument("-o", "--output", help="Directory for output")
    return parser.parse_args()


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
                    file_list.append([file_name, read_count])
                    bam_dict[sample]["file_list"].append([file_name, read_count])
                    bam_dict[sample]["average_read_count"] = sum(read_counts)/len(read_counts)
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


def return_subsample_data():
    bam_dict, file_list = generate_bam_dict()
    cut_list = calculate_clustering(bam_dict, 20)
    subsample_data = []

    for cut in cut_list:
        subsample = []
        for file in file_list:
            if file[1] <= cut:
                bam_file = [file[0], file[1], file[1]]
            else:
                bam_file = [file[0], file[1], cut]
            subsample.append(bam_file)
        subsample_data.append(subsample)

    return subsample_data


if __name__ == "__main__":

    CONFIG = read_config()
    args = parse_arguments()

    # Testing config
    bam_list = CONFIG['bam_list']

    # Parsing arguments

    # if args.cut is not None:
    #     cut = args.cut
    # else:
    #     test_folder = CONFIG['cut']
    #
    # if args.directory is not None:
    #     directory = args.directory
    # else:
    #     sys.stderr.write("Error: Must supply directory of bam files!")
    #     sys.exit()
    #
    # if args.output is not None:
    #     output = args.output
    # else:
    #     sys.stderr.write("Error: Must supply output directory!")
    #     sys.exit()
    #
    # if args.bam_list is not None:
    #     bam_list = args.bam_list
    # else:
    #     sys.stderr.write("Error: Must supply list of bam files and read counts!")
    #     sys.exit()

    data = return_subsample_data()
    sys.stdout.write(str(data))
