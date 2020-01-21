#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2019-10-03
# version: 	0.01

import sys
import re


# Main #
if __name__ == "__main__":

    usage = """
    ncbi-taxonomy.py <nodes.dmp> <names.dmp> <divison_ID>
    """

    if len(sys.argv) < 4:
        print(usage)
        sys.exit(1)

    nodesFile = sys.argv[1]
    namesFile = sys.argv[2]
    divId = sys.argv[3]
    namesObj = {}
    nodesObj = {}


    with open(namesFile) as names:

        for l in names:
            lA = re.split('\s\|\s', l.strip())
            lA[len(lA)-1] = lA[len(lA)-1].strip('\t|')
            if lA[0] not in namesObj:
                namesObj[lA[0]] = {
                    "ScientificName": "",
                    "SynonymName": "",
                    "CommonName": ""
                }
            # print(lA)
            if lA[-1] == "scientific name":
                namesObj[lA[0]]["ScientificName"] = lA[1]
            elif lA[-1] == "equivalent name":
                namesObj[lA[0]]["SynonymName"] += lA[1]
                namesObj[lA[0]]["SynonymName"] += ";"
            elif "synonym" in lA[-1]:
                namesObj[lA[0]]["SynonymName"] += lA[1]
                namesObj[lA[0]]["SynonymName"] += ";"
            elif "common" in lA[-1]:
                namesObj[lA[0]]["CommonName"] += lA[1]
                namesObj[lA[0]]["CommonName"] += ";"


    print("TaxID\tParentTaxId\tRank\tScientificName\tSynonymName\tCommonName")
    with open(nodesFile) as nodes:
        for l in nodes:

            lA = re.split('\s\|\s', l.strip())
            lA[len(lA)-1] = lA[len(lA)-1].strip('\t|')
            # print(lA)
            if lA[4] != divId:
                continue

            print("%s\t%s\t%s\t%s\t%s\t%s"%(lA[0],lA[1],lA[2],namesObj[lA[0]]["ScientificName"],namesObj[lA[0]]["SynonymName"],namesObj[lA[0]]["CommonName"]))
