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

    echo -ne "Downloading NCBI ARG Reference file $release [$outdir].."
    wget --quiet https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/Data/latest/ReferenceGeneCatalog.txt
    # python $script_path/ncbi_entrez.py bioproject_nucccore 313047 gb > NCBI-ARG.gbk
    echo " Done!"
}

# gunzip files
function organize_files() {
    echo -ne "Extracting CDS sequences from NCBI ..."
    script_path=`dirname $0`
    python $script_path/ncbi_arg.py ReferenceGeneCatalog.txt ncbi_arg.fasta
    echo " Done!"
}

download_files $1
organize_files
