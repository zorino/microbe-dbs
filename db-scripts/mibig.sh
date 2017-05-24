#!/bin/bash
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com


# download files from web site
function download_files() {

    release="1.3"

    out_dir=""
    if [ -z $1 ]
    then
        outdir=mibig_$release
    else
        outdir=$1/mibig_$release
    fi
    mkdir -p $outdir && cd $outdir

    echo "Downloading MIBiG $release [$outdir].."
    wget --quiet http://mibig.secondarymetabolites.org/mibig_json_$release.tar.gz
    wget --quiet http://mibig.secondarymetabolites.org/mibig_gbk_$release.tar.gz
    wget --quiet http://mibig.secondarymetabolites.org/MIBiG_prot_seqs_$release.fasta

}

# gunzip files
function organize_files() {
    for i in *.tar.gz
    do
        mkdir ${i%.tar.gz}
        tar xvf $i -C ${i%.tar.gz}
    done
    rm *.tar.gz
}


download_files $1
organize_files

