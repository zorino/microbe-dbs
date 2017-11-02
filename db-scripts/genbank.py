#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2017-05-25
# version: 	0.01

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def get_cds_seq(gb_file):

    with open(gb_file, "rU") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        header = record.id + "|"
                        header += feature.qualifiers['gene'][0] + "|"
                        header += feature.qualifiers['product'][0]
                        outseq = SeqRecord(Seq(feature.qualifiers['translation'][0],
                                               IUPAC.protein),
                                           id=header, name=feature.qualifiers['gene'][0],
                                           description="")
                        print(outseq.format("fasta").rstrip())



if __name__ == "__main__":

    usage = """
    genbank.py

         get_cds_seq           get cds sequences from genbank entries

    """

    if sys.argv[1] == "get_cds_seq":
        get_cds_seq(sys.argv[2])

    else:
        print(usage)
