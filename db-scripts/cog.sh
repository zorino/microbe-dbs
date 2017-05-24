#!/bin/bash
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com


# download files from web site
function download_files() {

    release=$(date +%Y-%m-%d)

    out_dir=""
    if [ -z $1 ]
    then
        outdir=cog_$release
    else
        outdir=$1/cog_$release
    fi
    mkdir -p $outdir && cd $outdir

    echo "Downloading COG $release [$outdir].."
    wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/prot2003-2014.fa.gz
    wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cog2003-2014.csv
    wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/prot2003-2014.gi2gbk.tab
    wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab
    wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab

}

# gunzip files
function organize_files() {
    gunzip *.gz
}


download_files $1
organize_files

