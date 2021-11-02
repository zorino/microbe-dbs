#!/bin/bash
# author : maxime déraspe
# email : maximilien1er@gmail.com
#
# required :
# - bash
# - wget

# download files from ftp
function download_ftp() {

    release="2019_09"

    out_dir=""
    if [ -z $1 ]
    then
        outdir=ebi_mgnify_uhgp90_$release
    else
        outdir=$1/ebi_mgnify_uhgp90_$release
    fi
    mkdir -p $outdir && cd $outdir

    echo "Downloading EBI Mgnify Protein Catalogue $release [$outdir].."
    wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgp_catalogue/uhgp-90.tar.gz

}

# organize files into sub-directories
function organize_files () {

    echo -e "\n\nMerging Mgnify Data.."
    dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    tar -zxvf uhgp-90.tar.gz uhgp-90/uhgp-90_hq.faa uhgp-90/uhgp-90_eggNOG.tsv
    python $dir/ebi_mgnify_uhgp.py uhgp-90/uhgp-90_hq.faa uhgp-90/uhgp-90_eggNOG.tsv > uhgp-90.tsv
    rm -fr uhgp-90.tar.gz uhgp-90/

}


download_ftp $1
organize_files
