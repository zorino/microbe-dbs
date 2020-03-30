#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2020-03-27
# version: 	0.01

import sys
import csv


def read_fasta(fasta_file, csv_header, csv):

    print("\t".join(csv_header))
    sequences = {}
    with open(fasta_file) as f:
        for l in f:
            if l[0] == ">":
                lA = l.split()
                id = lA[0].replace(">", "")
                orf = l.split("~~~")[2]
                sequences[id] = {"gene": orf, "seq": ""}
            else:
                sequences[id]["seq"] += l.rstrip()

    for k, v in sequences.items():
        csv[k][1] += " %s" % v["gene"]
        csv[k].append(v["seq"])
        print("\t".join(csv[k]))


def read_csv(csv_file):

    cols_to_keep = [0, 1, 2, 4, 6, 7, 8, 9, 10, 11, 12]
    f = open(csv_file, 'r')

    csv_table = {}

    header = f.readline().rstrip().split(",")
    header_sub = [header[i] for i in cols_to_keep]
    header_sub[0] = "EntryID"
    header_sub[1] = "ProteinName"
    header_sub.append("Sequence")

    reader = csv.reader(f, delimiter=',')
    i = 0
    for row in reader:
        row[0] = "%s" % i
        output = [row[i] for i in cols_to_keep]
        # print('\t'.join(output))
        i += 1
        csv_table[output[1]] = output

    return header_sub, csv_table


# Main #
if __name__ == "__main__":
    usage = """
    python isfinder.py <fasta> <csv>
"""

    if len(sys.argv) < 3:
        print(usage)
        sys.exit(1)

    csv_header, csv = read_csv(sys.argv[2])
    sequences = read_fasta(sys.argv[1], csv_header, csv)
