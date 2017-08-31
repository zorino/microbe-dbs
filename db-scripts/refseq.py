#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2017-08-24
# version: 	0.01

from subprocess import call
import argparse
import gzip
import re

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


def load_taxonomy(nodes_file, division):

    taxa = {}
    with open(nodes_file) as n:
        for l in n:
            fields = re.compile("[\s+\t+]\|[\s+\t+]").split(l.strip())
            if fields[4] != division:
                continue
            if fields[1] not in taxa:
                taxa[fields[1]] = {'child': [], 'parent': "", 'rank': "", 'path': ""}
            if fields[0] not in taxa:
                taxa[fields[0]] = {'child': [], 'parent': "", 'rank': "", 'path': ""}
            taxa[fields[1]]['child'].append(fields[0])
            taxa[fields[0]]['parent'] = fields[1]
            taxa[fields[0]]['rank'] = fields[2]
            # print(fields)

    return taxa


def create_ancestor_path(taxa, tax_id):

    path = []
    path.append(tax_id)
    while taxa[tax_id]['parent'] != "":
        tax_id = taxa[tax_id]['parent']
        path.insert(0, tax_id)

    return path

def create_taxonomy_directories(taxa):

    taxon_done = []

    for t in taxa:

        path = create_ancestor_path(taxa, t)
        path_str = "/".join(path)
        taxa[t]['path'] = path_str

        if t not in taxon_done:
            taxon_done.extend(path)
            print("create dir %s" % path_str)
            call(['mkdir', '-p', path_str])

        with open(path_str+"/rank", "w") as rank_file:
            rank_file.write(taxa[t]['rank'])

    return 0

def create_names(taxa, names_file):

    with open(names_file) as n:
        for l in n:
            fields = re.compile("[\s+\t+]\|[\s+\t+]").split(l.strip())
            if fields[3] == "scientific name":
                with open(taxa[fields[0]]['path']+"/name", "w") as name_file:
                    name_file.write(fields[1])
            else:
                with open(taxa[fields[0]]['path']+"/name_other", "w") as name_file:
                    out = fields[1] + "\t" + fields[2]
                    name_file.write(out)

    return 0

def taxonomy(opt):

    taxa = load_taxonomy(opt.nodes, opt.div)
    create_taxonomy_directories(taxa)
    if opt.names:
        create_names(taxa, opt.names)

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

    # subparser to search sequence
    parser_search = subparsers.add_parser('taxonomy', help='group sequence by taxonomic id inside directories')
    parser_search.add_argument("--ids", type=str, required=True, metavar="IDS-TAXA")
    parser_search.add_argument("--nodes", "-n", type=str, required=True, metavar="nodes.dmp")
    parser_search.add_argument("--names", type=str, required=True, metavar="names.dmp")
    parser_search.add_argument("--div", "-d", type=str, required=True, metavar="bacteria, viruses, ..")
    parser_search.add_argument("--fna", "-f", type=str, required=False)
    parser_search.set_defaults(which='taxonomy')

    args = parser.parse_args()

    if args.which == "search":
        search_feature(args)
    elif args.which == "taxonomy":
        taxonomy(args)
