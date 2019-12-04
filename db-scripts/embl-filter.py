#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2019-09-03
# version: 	0.01

import sys
import gzip
import re

# Get only taxon entries from EMBL input file
def taxon(taxon, embl_file):

    if embl_file[-2:] == "gz":
        filein = gzip.open(embl_file, "rt")
    else:
        filein = open(embl_file)

    entry = ""
    keep = 0

    for l in filein:
        if l[0:2] == "//":
            if keep == 1:
                entry += l
                print(entry)
            entry = ""
            keep = 0
        else:
            entry += l

            if l[0:2] == "OS" or l[0:2] == "OC":
                if taxon in l:
                    keep = 1

    return

def fasta(embl_file):

    if embl_file[-2:] == "gz":
        filein = gzip.open(embl_file, "rt")
    else:
        filein = open(embl_file)

    entry = ""
    keep = 0
    reg = re.compile("\s+")
    inside_seq = False

    for l in filein:
        if l[0:2] == "//":
            if entry != "":
                print(entry)
            entry = ""
            inside_seq = False
        elif l[0:2] == "ID":
            entry += ">%s\n" % (reg.split(l))[1]
        elif l[0:2] == "SQ":
            inside_seq = True
        elif inside_seq:
            seq_split = reg.split(l.strip())
            entry += "".join(seq_split)

    return


# Main #
if __name__ == "__main__":

    usage = """

embl-filter.py [options] INPUT_EMBL

  taxon          <taxon> (ex. Escherichia coli)

  fasta          extract fasta sequence


    """

    if len(sys.argv) < 2:
        print(usage)
        exit(1)

    if sys.argv[1] == "taxon":
        if len(sys.argv) < 4:
            print(usage)
            exit(1)

        taxon(sys.argv[2], sys.argv[-1])

    elif sys.argv[1] == "fasta":
        if len(sys.argv) < 3:
            print(usage)
            exit(1)
        fasta(sys.argv[2])

    else:
        print(usage)
        exit(1)
