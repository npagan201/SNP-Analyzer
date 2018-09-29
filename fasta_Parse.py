#!/usr/bin/env python3

import sys


def fastaParse(infile):
    with open(infile, 'r') as fastaFile:
        # Skip whitespace
        while True:
            line = fastaFile.readline()
            if line is "":
                return  # Empty file or premature end of file?
            if line[0] is ">":
                break
        while True:
            if line[0] is not ">":
                raise ValueError("Records in FASTA should begin with '>'")
            header = line[1:].rstrip()
            allLines = []
            line = fastaFile.readline()
            while True:
                if not line:
                    break
                if line[0] is ">":
                    break
                allLines.append(line.rstrip())
                line = fastaFile.readline()
            yield header, "".join(allLines).replace(" ", "").replace("\r", "")
            if not line:
                return  # Stop Iteration


def load_annotation_file(infile):
    return_values = {}
    with open(infile, 'r') as f:
        data = f.read().split('\n')
        for line in data:
            if not line:
                continue
            entries = line.split(',')
            return_values.setdefault(entries[0], entries[1])
    return return_values


def output_data(D, A, outfile):
    with open(outfile, 'w') as aout:
        aout.write('Header,Annotations\n')
        for header, seq in D.items():
            if header in A:
                sys.stdout.write('>{}\n{}\n'.format(
                    header.replace(' ', '_'),
                    seq
                ))
                aout.write('{},\n'.format(
                    header.replace(' ', '_')
                ))


if __name__ == '__main__':
    D = {k: v for k, v in fastaParse(sys.argv[1])}
    A = load_annotation_file(sys.argv[2])
    output_data(D, A, sys.argv[3])