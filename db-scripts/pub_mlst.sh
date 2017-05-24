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
        outdir=pubmlst_$release
    else
        outdir=$1/pubmlst_$release
    fi
    mkdir -p $outdir && cd $outdir

    echo "Downloading PubMLST $release [$outdir].."
    wget --quiet https://pubmlst.org/data/dbases.xml

}

# gunzip files
function organize_files() {

    echo "Parsing / Loading Locus from PubMLST $release.."
    script_path=`dirname $0`
    python $script_path/pub_mlst.py dbases.xml > dbases.tsv

    # download all mlst locus files
    IFS=$'\t'
    while read a b c d e
    do
        species=`echo $c | sed "s#https://pubmlst.org/data/profiles/##g"`
        echo ${species%.txt} $d
        if [ ! -f $species ]
        then
            wget --quiet $c
            mkdir ${species%.txt}
        fi
        cd ${species%.txt}
        wget --quiet $e
        rename .tfa .fasta *
        cd ../
    done < dbases.tsv

}


download_files $1
organize_files

