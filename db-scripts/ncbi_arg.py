#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2019-10-09
# version: 	0.01

import sys
import time
from Bio import Entrez

def ncbi_fetch_protein(protein_id):

    Entrez.email = "microbesdbs@example.com"
    handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
    data = handle.read()
    print(data.strip())


# Main #
if __name__ == "__main__":

    file = sys.argv[1]

    with open(file) as f:
        l = f.readline()
        for l in f:
            lA = l.split("\t")
            if lA[9] != "":
                ncbi_fetch_protein(lA[9])
            elif lA[12] != "":
                ncbi_fetch_protein(lA[12])

            time.sleep(2)
