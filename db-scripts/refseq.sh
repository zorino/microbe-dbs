#!/bin/bash
# author : maxime déraspe
# email : maximilien1er@gmail.com
#
# required :
# - bash
# - wget

# download files from ftp
function download_ftp() {

    release=`wget -O- -q ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER`;

    out_dir=""
    if [ -z $1 ]
    then
        outdir=refseq_$2_$release
    else
        outdir=$1/refseq_$2_$release
    fi
    mkdir -p $outdir && cd $outdir

    echo "Downloading NCBI RefSeq $2 $release [$outdir].."
    wget -R --no-parent  ftp://ftp.ncbi.nlm.nih.gov/refseq/release/$2/*.genomic.*
    wget -R --no-parent  ftp://ftp.ncbi.nlm.nih.gov/refseq/release/$2/*.rna.*
    wget -R --no-parent  ftp://ftp.ncbi.nlm.nih.gov/refseq/release/$2/*protein.*

}

# organize files into sub-directories
function organize_files () {

    files=(*protein.*)
    if [ -e "${files[0]}" ]
    then
        mkdir protein
        mv *protein.* protein
    fi

    files=(*.genomic.*)
    if [ -e "${files[0]}" ]
    then
        mkdir genomic
        mv *.genomic.* genomic
    fi

    files=(*.rna.*)
    if [ -e "${files[0]}" ]
    then
        mkdir rna
        mv *.rna.* rna
    fi

}


download_ftp $1 $2
organize_files

