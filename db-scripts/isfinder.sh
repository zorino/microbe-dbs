#!/bin/bash
# author : maxime déraspe
# email : maximilien1er@gmail.com
#
# required :
# - bash
# - wget

# download files from ftp
function download_ftp() {

    release=2019-09-25

    out_dir=""
    if [ -z $1 ]
    then
        outdir=isfinder_$release
    else
        outdir=$1/isfinder_$release
    fi
    mkdir -p $outdir && cd $outdir

    echo "Downloading Isfinder $2 $release [$outdir].."
    wget https://github.com/zorino/ISfinder-sequences/archive/$release.tar.gz -O ISfinder.tgz

}

# organize files into sub-directories
function organize_files () {
    tar zxf ISfinder.tgz
    rm ISfinder.tgz
    mv ISfinder-sequences-2019-09-25/IS* .
    rm -fr ISfinder-sequences-2019-09-25
    dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    python $dir/isfinder.py IS.faa IS.csv > IS.tsv
}


download_ftp $1 $2
organize_files
