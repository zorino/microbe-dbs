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
        outdir=vfdb_$release
    else
        outdir=$1/vfdb_$release
    fi
    mkdir -p $outdir && cd $outdir

    echo "Downloading VFDB $release [$outdir].."
    wget --quiet http://www.mgc.ac.cn/VFs/Down/VFs.xls.gz
    wget --quiet http://www.mgc.ac.cn/VFs/Down/Comparative_tables_from_VFDB.tar.gz
    wget --quiet http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz
    wget --quiet http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
    wget --quiet http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz
    wget --quiet http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz

}

# gunzip files
function organize_files() {
    gunzip *.gz
    tar xvf *.tar
    rm *.tar
}


download_files $1
organize_files

