#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2019-10-09
# version: 	0.01

import sys
import time
from Bio import Entrez

def ncbi_fetch_protein(protein_id, output):

    Entrez.email = "microbesdbs@example.com"
    handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
    data = handle.read()
    output.write(data.strip())
    output.write("\n")


# Main #
if __name__ == "__main__":

    file = sys.argv[1]
    file_output = sys.argv[2]

    output = open(file_output, "w")

    counter = 0

    with open(file) as f:
        l = f.readline()
        for l in f:
            counter += 1
            lA = l.split("\t")
            if lA[9] != "":
                ncbi_fetch_protein(lA[9], output)
            elif lA[12] != "":
                ncbi_fetch_protein(lA[12], output)

            if (counter % 2) == 0:
                time.sleep(1)


    output.close
