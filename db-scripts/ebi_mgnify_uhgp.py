#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2020-02-12
# version: 	0.01

import sys


def merge_file(fasta_file, eggnog_file):

    fasta = {}
    seq = ""
    id = ""
    name = ""

    with open(fasta_file) as f:
        for l in f:

            if l[0] == ">":
                if seq != "":
                    fasta[id] = {"name": name, "seq": seq}

                id = l.strip().split(" ")[0].replace(">", "")
                name = " ".join(l.strip().split(" ")[1:])
                seq = ""
            else:
                seq += l.strip()

        if seq != "":
            fasta[id] = {"name": name, "seq": seq}

    print(
        "EntryID\tProteinName\tGeneName\tTaxon.Class\tGO\tEC\tKEGG.orth\tKEGG.pathway\tKEGG.module\tCOG.function\tSequence"
    )

    with open(eggnog_file) as f:
        for l in f:
            lA = l.rstrip('\n').split("\t")
            kegg_pathways = ""
            for e in lA[9].split(","):
                if kegg_pathways != "":
                    kegg_pathways += (",%s" % e)
                else:
                    kegg_pathways += ("%s" % e)
            cog_function = ""
            if len(lA) > 20:
                cog_function = lA[20]
            if lA[0] in fasta:
                print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                      (lA[0], fasta[lA[0]]["name"], lA[5], lA[4], lA[6], lA[7],
                       lA[8], kegg_pathways, lA[10], cog_function,
                       fasta[lA[0]]["seq"]))


# Main #
if __name__ == "__main__":

    usage = """
    ebi-mgnifyh-uhgp.py <fasta> <tsv>
"""

    if len(sys.argv) < 3:
        print(usage)
        sys.exit(1)

    merge_file(sys.argv[1], sys.argv[2])
