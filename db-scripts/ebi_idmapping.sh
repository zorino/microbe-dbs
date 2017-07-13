#!/bin/bash
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com


# download files from web site
function download_files() {

    # release=$(date +%Y-%m-%d)
    release=`wget -O- -q ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/RELEASE.metalink | grep \<version | cut -d ">" -f2 | cut -d "<" -f1`;

    out_dir=""
    if [ -z $1 ]
    then
        outdir=ebi-idmapping_$release
    else
        outdir=$1/ebi-idmapping_$release
    fi
    mkdir -p $outdir && cd $outdir

    echo "Downloading VFDB $release [$outdir].."
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/RELEASE.metalink

}

# gunzip files
function organize_files() {
    echo "organizing files"
}


download_files $1
# organize_files

