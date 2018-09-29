#!/usr/bin/env python3

import os
import sys


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
                    return temp[1], temp[2], temp[3], temp[9]

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
    set()
    values = S.next()
    while values:
        print(values)
        try:
            values = S.next()
        except StopIteration:
            break
