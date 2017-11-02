#!/bin/bash
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com


# download files from web site
function download_files() {

    release=$(date +%Y-%m-%d)
	script_path=`dirname $0`

    outdir=""
    if [ -z $1 ]
    then
        outdir=ncbi-arg_$release
    else
        outdir=$1/ncbi-arg_$release
    fi
    mkdir -p $outdir && cd $outdir

    echo -ne "Downloading NCBI ARG $release [$outdir].."
	python $script_path/ncbi_entrez.py 313047 gb > NCBI-ARG.gbk
	echo " Done!"
}

# gunzip files
function organize_files() {
	echo -ne "Extracting CDS sequences from genbank.."
	script_path=`dirname $0`
	python $script_path/genbank.py get_cds_seq NCBI-ARG.gbk > NCBI-ARG.fasta
	echo " Done!"
}

download_files $1
organize_files
