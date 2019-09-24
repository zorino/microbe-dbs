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
        outdir=uniprotkb_bacteria_$release
    else
        outdir=$1/uniprotkb_bacteria_$release
    fi
    mkdir -p $outdir && cd $outdir

    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_bacteria.dat.gz
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_bacteria.dat.gz

}


function organize_files() {

    echo "Organizing files - creating fasta"
    script_path=`dirname $0`
    python $script_path/embl-filter.py fasta uniprot_sprot_bacteria.dat.gz > uniprot_sprot_bacteria.fasta
    python $script_path/embl-filter.py fasta uniprot_trembl_bacteria.dat.gz > uniprot_trembl_bacteria.fasta

}


download_files $1
organize_files
