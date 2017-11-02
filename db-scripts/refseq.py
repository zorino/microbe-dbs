#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2017-08-24
# version: 	0.01

from subprocess import call
import os
import argparse
import gzip
import re
import glob
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


def create_ancestor_path(taxa, tax_id):

    path = []
    path.append(tax_id)
    while taxa[tax_id]['parent'] != "":
        tax_id = taxa[tax_id]['parent']
        path.insert(0, tax_id)

    return path

def create_ancestor_path_all(taxa):

    for t in taxa:
        path = create_ancestor_path(taxa, t)
        path_str = "/".join(path)
        taxa[t]['path'] = path_str

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

    create_ancestor_path_all(taxa)

    return taxa


def create_taxonomy_directories(taxa):

    taxon_done = []

    for t in taxa:

        path_str = taxa[t]['path']
        path = path_str.split("/")

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

            fields = re.compile("[\s*\t*]\|[\s*\t*]").split(l)

            if fields[0] not in taxa:
                continue

            print("writing names into %s" % taxa[fields[0]]['path'])

            if fields[3] == "scientific name":
                with open(taxa[fields[0]]['path']+"/name", "w") as name_file:
                    name_file.write(fields[1])
            else:
                with open(taxa[fields[0]]['path']+"/name_other", "a") as name_file:
                    out = fields[1] + "\t" + fields[3] + "\n"
                    name_file.write(out)

    return 0

def load_ids(ids_file, taxa):

    genomes = {}
    ids_files = ids_file.split(",")
    for ids_f in ids_files:
        print(ids_f)
        ids = open_file(ids_f)
        for l in ids:
            fields = re.compile("\t").split(l.strip())
            if fields[2] in taxa:
                genomes[fields[0]] = fields[2]
        ids.close()

    return genomes

def create_sequences(genome_ids, taxa, fna_dir, gbf_dir):

    all_fna = glob.glob(fna_dir + "/*.fna.gz")

    all_gbf = glob.glob(gbf_dir + "/*.gbff.gz")

    for gbf in all_gbf:
        with open_file(gbf) as f:
            rec_gbf = SeqIO.parse(f, "genbank")
            file_id = ".".join(gbf.split("/")[-1:][0].split(".")[:2]) + "."
            file_fna = [s for s in all_fna if file_id in s][0]
            handle_fna = open_file(file_fna)
            rec_fna = SeqIO.parse(handle_fna, "fasta")
            for r_gbf in rec_gbf:

                r_fna = next(rec_fna)
                biosample = ""
                biosample_xref = [s for s in r_gbf.dbxrefs if "BioSample" in s]
                biosample_assembly = [s for s in r_gbf.dbxrefs if "Assembly" in s]
                if len(biosample_xref) > 0:
                    biosample= biosample_xref[0].replace("BioSample:","")
                elif len(biosample_assembly) > 0:
                    biosample = biosample_assembly[0].replace("Assembly:","")
                else:
                    biosample = "OTHER_SAMPLES"
                    # continue

                print("%s\t%s\t%s" % (biosample, r_gbf.id, ";".join(r_gbf.dbxrefs)))

                genome_id = r_gbf.id.split(".")[0]
                if genome_id in genome_ids:
                    path_str = taxa[genome_ids[genome_id]]['path']
                else:
                    path_str = "./unidentified"

                samples_path = path_str + "/samples"

                if not os.path.exists(samples_path):
                    call(['mkdir', '-p', samples_path ])

                print(samples_path)

                with open(samples_path+"/"+biosample+".fasta","a") as f:
                    f.write(r_fna.format("fasta"))

                with open(samples_path+"/"+biosample+".gb","a") as f:
                    f.write(r_gbf.format("genbank"))

            handle_fna.close()

    return 0


def taxonomy(opt):

    print("Loading taxonomy nodes..")
    taxa = load_taxonomy(opt.nodes, opt.div)

    print("Creating taxonomy directories..")
    create_taxonomy_directories(taxa)

    if opt.names:
        print("Creating names for each taxon..")
        create_names(taxa, opt.names)

    print("Loading genome ids..")
    genome_ids = load_ids(opt.ids, taxa)

    print("Dumping sequences to taxon directory..")
    create_sequences(genome_ids, taxa, opt.fna, opt.fna)

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
    parser_search = subparsers.add_parser('taxonomy', help='group sequence by taxonomic ids inside directories')
    parser_search.add_argument("--ids", type=str, required=True, metavar="IDS-TAXA")
    parser_search.add_argument("--nodes", "-n", type=str, required=True, metavar="nodes.dmp")
    parser_search.add_argument("--names", type=str, metavar="names.dmp")
    parser_search.add_argument("--div", "-d", type=str, required=True, metavar="bacteria(0), phages(3), viruses(9)..")
    parser_search.add_argument("--fna", "-f", type=str, required=False, metavar="FNA_Directory")
    parser_search.set_defaults(which='taxonomy')

    args = parser.parse_args()

    if hasattr(args, "which"):
        if args.which == "search":
            search_feature(args)
        elif args.which == "taxonomy":
            taxonomy(args)
        else:
            parser.print_help()

    else:
        # print(ValueError)
        parser.print_help()
