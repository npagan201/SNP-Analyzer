#!/usr/bin/env python3

import csv
import os
import sys
import re

"""This method pulls out the headers from a given AMR_matrix and a SNP_metadata file and 
pulls out the common headers between the two"""


def header_collect(amr_matrix, snp_metadata):
    matrix_headers = {}
    snp_headers = {}
    temp_list = []
    with open(amr_matrix, 'r') as csvfile:  # pulls the headers from the AMR_matrix, delimiting by comma
        reader = csv.reader(csvfile, delimiter=',')
        csvfile.readline()
        for row in reader:
            if row != []:
                matrix_headers[row[0]] = ''
    with open(snp_metadata, 'r') as csvfile: # pulls the headers from the SNP_metadata, delimiting by comma
        reader = csv.reader(csvfile, delimiter=',')
        csvfile.readline()
        for row in reader:
            if row != []:
                if matrix_headers.__contains__(row[0]):  # pulls out the headers if they are contained in the AMR_matrix
                    temp_list.append([row[2], row[3], row[4], row[5]]) # row[10]] # pulls out start/stop index & wild/mutant type amino acid & sequence
                    snp_headers[row[0]] = temp_list
    return snp_headers


def cigar_str_to_int(cig_list):
    """Convert int strings in the cigar_list to ints"""
    i = 0
    for x in cig_list[::2]:
        cig_list[i] = int(x)
        i += 2
    return cig_list


def cigar_count(cig_list):
    """Count of how long the sequence should be based on the given cigar string"""
    i = 0
    for x in cig_list[::2]:
        i += x
    return i


def start_stop_to_int(dictionary):
    """Convert start/stop indices given by the SNP metadata to int instead of str,
    if the indices of the values are 0 and 1"""
    for lists in dictionary:
        for each_list in range(len(dictionary[lists])):
            for start_stop in range(0, 2):
                dictionary[lists][each_list][start_stop] = int(dictionary[lists][each_list][start_stop])
    return dictionary


def cigar_switch(operator):
    switch = {
        'M': "M_here",
        'X': "X_here",
        '=': "=_here",
        'D': "D_here",
        'N': "N_here",
        'S': "S_here",
        'I': "I_here",
        'H': "H_here",
        'P': "P_here",
    }
    print(switch.get(operator, "Invalid cigar operator"))


class SamParser:
    """This object takes as input a SAM file path and constructs an iterable that outputs
     sequence information.  Only one line will be held in memory at a time using this method.
    """
    def __init__(self, filepath):
        """
        constructor
        @param filepath: filepath to the input raw SAM file.
        """
        if os.path.exists(filepath):  # if file is a file, read from the file
            self.sam_file = str(filepath)
            self.stdin = False
        elif not sys.stdin.isatty():  # else read from standard in
            self.stdin = True
        else:
            raise ValueError("Parameter filepath must be a SAM file")
        self.current_line = None
        self.reads_mapping = 0
        self.reads_total = 0
        # Allowable bitflags for SAM file -> reads with both mates mapping, regardless of other flags
        self.true_flags = (99, 147, 83, 163, 67, 131, 115, 179, 81, 161, 97, 145, 65, 129, 113, 177)

    def __iter__(self):
        return self

    def _iterate(self):
        # Skip all leading whitespace
        while True:
            if self.stdin:
                sam_line = sys.stdin.readline()  # read from stdin
            else:
                sam_line = self.sam_file.readline()  # read from file
            if not sam_line:
                return  # End of file
            if sam_line.startswith("@SQ"):  # these lines contain refseq length information
                temp = sam_line.split()
                return temp[1][3:], temp[2][3:]
            elif sam_line[0] != "@":  # these lines are the actual reads
                self.reads_total += 1
                if self.reads_total % 100000 == 0:  # update the counter on stdout every 100000 reads
                    sys.stdout.write("\rReads processed: {}".format(self.reads_total))
                    sys.stdout.flush()
                temp = sam_line.split()
                if int(temp[1]) in self.true_flags and temp[2] is not "*" and int(temp[3]) is not 0:
                    self.reads_mapping += 1
                    return temp[1], temp[2], temp[3], temp[5], temp[9]

    def next(self):
        if not self.stdin and type(self.sam_file) is str:  # only open file here if sam_file is a str and not file
            self.sam_file = open(self.sam_file, "r")
        value = self._iterate()
        if not value:  # close file on EOF
            if not self.stdin:
                self.sam_file.close()
            print("{0} with both mates mapped out of {1} total reads\n".format(self.reads_mapping, self.reads_total))
            raise StopIteration()
        else:
            return value


if __name__ == '__main__':
    S = SamParser(sys.argv[1])
    snp_data = header_collect(sys.argv[2], sys.argv[3])  # call the method above with command line arguments
    snp_data = start_stop_to_int(snp_data)
    set()
    values = S.next()
    while values:
        if snp_data.__contains__(values[1]):  # checks if the dictionary contains the given header
            cigar_list = list(filter(None, re.split('(\d+)', values[3])))  # get cigar string and divide it into numbers/letters
            cigar_list = cigar_str_to_int(cigar_list)
            cigar_total = cigar_count(cigar_list)
            for x in snp_data:
                for y in range(0, len(snp_data[x])):
                    if (int(values[2]) >= snp_data[x][y][0]) and (int(values[2]) <= snp_data[x][y][1]):
                        print(values[2])
        try:
            values = S.next()
        except StopIteration:
            break

