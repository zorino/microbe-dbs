#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2020-03-30
# version: 	0.01

import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def read_fasta(files):

    sequences = {}
    for file in files:
        with open(file) as f:
            for l in f:
                if l[0] == ">":
                    lA = l.split()
                    id = lA[0].replace(">", "")
                    sequences[id] = ""
                else:
                    sequences[id] += l.rstrip()

    return sequences


def parse_pheno(pheno_file, sequences):

    with open(pheno_file) as f:
        _header_str = f.readline()
        header = _header_str.rstrip().split("\t")
        header[0] = "EntryID"
        header.append("ProteinName")
        header.append("Sequence")
        print("\t".join(header))

        for l in f:
            lA = l.rstrip("\n").split("\t")
            if lA[0] not in sequences:
                continue
            protein_name = lA[0].split("_")
            protein_name = "_".join(protein_name[:len(protein_name) - 1])
            lA.append(protein_name)
            lA.append(
                str(
                    Seq(sequences[lA[0].strip()],
                        generic_dna).translate(table=11)))
            print("\t".join(lA))


# Main #
if __name__ == "__main__":

    usage = """
    python resfinder_db.py <phenotypes.txt> *.fsa
"""

    if len(sys.argv) < 3:
        print(usage)
        sys.exit(1)

    sequences = read_fasta(sys.argv[2:])

    parse_pheno(sys.argv[1], sequences)
