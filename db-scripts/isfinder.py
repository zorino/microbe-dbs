#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2020-03-27
# version: 	0.01

import sys
import csv


def read_fasta(fasta_file):

    sequences = {}
    with open(fasta_file) as f:
        for l in f:
            if l[0] == ">":
                lA = l.split()
                id = lA[0].replace(">", "")
                sequences[id] = ""
            else:
                sequences[id] += l.rstrip()

    return sequences


def read_csv(csv_file, sequences):

    cols_to_keep = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13]
    f = open(csv_file, 'r')

    header = f.readline().rstrip().split(",")
    header.append("Sequence")
    header_sub = [header[i] for i in cols_to_keep]
    print("\t".join(header_sub))

    reader = csv.reader(f, delimiter=',')
    for row in reader:
        if row[1] in sequences:
            seq = sequences[row[1]]
            row.append(seq)
            output = [row[i] for i in cols_to_keep]
            print('\t'.join(output))


# Main #
if __name__ == "__main__":
    usage = """
    python isfinder.py <fasta> <csv>
"""

    if len(sys.argv) < 3:
        print(usage)
        sys.exit(1)

    sequences = read_fasta(sys.argv[1])
    csv = read_csv(sys.argv[2], sequences)

    # f = open(sys.argv[1])

    # reader = csv.reader(f.readlines(), delimiter=',')
    # for row in reader:
    #     print('\t'.join(row))
