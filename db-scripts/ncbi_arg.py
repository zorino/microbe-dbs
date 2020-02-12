#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2019-10-09
# version: 	0.01

import sys
import time
from Bio import Entrez, SeqIO
import urllib

Entrez.email = "microbesdbs@example.com"

def ncbi_fetch_proteins(protein_ids, output):

    request = Entrez.epost("protein", id=",".join(protein_ids))

    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        print("An error occurred while retrieving the annotations.")
        print("The error returned was %s" % e)
        sys.exit(-1)

    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.efetch(db="protein", webenv=webEnv, query_key=queryKey, retmode="xml")
    annotations = Entrez.read(data)

    for a in annotations:
        output.write(">%s\n"%a['GBSeq_accession-version'])
        output.write("%s\n"%a['GBSeq_sequence'].upper())


def merge_reference_with_sequence(tsv_file, fasta_file):

    new_output = tsv_file.replace(".txt",".tsv")
    writer = open(new_output, "w")

    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)

    with open(tsv_file) as f:
        header = f.readline().rstrip().split("\t")
        del header[14:19]
        del header[11]
        del header[2]
        header.insert(0,"entryid")
        header.append("sequence")
        header.append("sequence_length")
        writer.write("\t".join(header))
        writer.write("\n")

        for l in f:
            lA = l.rstrip().split("\t")
            entryid = ""
            if lA[9] != "":
                entryid = lA[9]
            elif lA[12] != "":
                entryid = lA[12]

            if entryid == "" or (entryid not in sequences):
                continue

            # delete unwanted index
            del lA[14:19]
            del lA[11]
            del lA[2]
            lA.insert(0, entryid)
            lA.append(sequences[entryid])
            lA.append(str(len(sequences[entryid])))
            writer.write("\t".join(lA))
            writer.write("\n")

# Main #
if __name__ == "__main__":

    file_meta = sys.argv[1]
    file_output = sys.argv[2]

    output = open(file_output, "w")

    counter = 0
    protein_ids = []

    with open(file_meta) as f:
        l = f.readline()
        for l in f:
            counter += 1
            lA = l.split("\t")
            if lA[5] != "AMR" and lA[5] != "STRESS":
                continue
            if lA[9] != "":
                if len(protein_ids) < 100:
                    protein_ids.append(lA[9])
            elif lA[12] != "":
                if len(protein_ids) < 100:
                    protein_ids.append(lA[12])

            if len(protein_ids) >= 100:
                ncbi_fetch_proteins(protein_ids, output)
                protein_ids = []


    if len(protein_ids) > 0:
        ncbi_fetch_proteins(protein_ids, output)

    output.close

    merge_reference_with_sequence(file_meta, file_output)
