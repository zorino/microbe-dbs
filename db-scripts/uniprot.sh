#!/bin/bash
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com


function download_files() {

	echo "Download sequences"

	release=$(date +%Y-%m-%d)

    out_dir=""
    if [ -z $1 ]
    then
        outdir=vfdb_$release
    else
        outdir=$1/vfdb_$release
    fi
    mkdir -p $outdir && cd $outdir

	# tsv file
	curl 'https://www.uniprot.org/uniprot/?query=taxonomy:2&format=tab&force=true&columns=id,entry%20name,reviewed,protein%20names,genes,organism,lineage(all),go,comment(FUNCTION),comment(PATHWAY),ec,mass,length,sequence&sort=score&compress=yes' > uniprotkb-bacteria-metadata.tsv.gz

}


function organize_files() {

	zcat uniprotkb-bacteria-metadata.tsv.gz  | awk -F "\t" '{print ">"$1"|"$4"\n"$14}' > uniprotkb-bacteria-sequence.fasta

}


download_files $1
organize_files
