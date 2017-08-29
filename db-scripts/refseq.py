#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2017-08-24
# version: 	0.01

import argparse
import gzip
from Bio import SeqIO

#from Bio.SeqFeature import FeatureLocation


def open_file(f):

    if f[-3:] == ".gz":
        handle = gzip.open(f, "rt")
    else:
        handle = open(f, "rt")

    return handle


# search for a feature
def search_feature(opt):

    handle_gbf = open_file(opt.genbank)
    handle_fna = open_file(opt.fna)

    rec_gbf = SeqIO.parse(handle_gbf, "genbank")
    rec_fna = SeqIO.parse(handle_fna, "fasta")

    for r_gbf in rec_gbf:
        r_fna = next(rec_fna)
        for ft in r_gbf.features:
            if opt.ft_type in ft.type:
                if opt.complete and (">" in str(ft.location) or "<" in str(ft.location)):
                    continue
                seq = r_fna.seq[ft.location.start:ft.location.end]
                if ft.location.strand != 1:
                    seq = seq.complement()
                print(">%s %s" %(r_gbf.id, ft.qualifiers['product'][0]))
                print(seq)

    handle_gbf.close()
    handle_fna.close()

    return 0


# Main #
if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='refseq.py')
    subparsers = parser.add_subparsers()

    # subparser to search sequence
    parser_search = subparsers.add_parser('ft_search', help='search for a feature in genbank file')
    parser_search.add_argument("--ft-type", "-t", type=str, default="", metavar="FEATURE_TYPE")
    parser_search.add_argument("--genbank", "-g", type=str, required=True)
    parser_search.add_argument("--fna", "-f", type=str, required=True)
    parser_search.add_argument("--complete", action="store_true", help="need a complete sequence, default: false")
    parser_search.set_defaults(which='search')


    args = parser.parse_args()

    if args.which == "search":
        search_feature(args)
