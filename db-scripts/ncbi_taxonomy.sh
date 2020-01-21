#!/bin/bash
# author : maxime déraspe
# email : maximilien1er@gmail.com
#
# required :
# - bash
# - wget

# download files from ftp
function download_ftp() {

    release=$(date +%Y-%m-%d)

    out_dir=""
    if [ -z $1 ]
    then
        outdir=ncbi-taxonomy_$release
    else
        outdir=$1/ncbi-taxonomy_$release
    fi
    mkdir -p $outdir && cd $outdir

    echo "Downloading NCBI Taxonomy $release [$outdir].."
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz

}

# organize files into sub-directories
function organize_files () {

    mkdir taxdump
    tar xf taxdump.tar.gz -C ./taxdump
    echo -e "\n\nExtracting Bacteria Kingdom Taxonomy.."
    dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    python $dir/ncbi_taxonomy.py ./taxdump/nodes.dmp ./taxdump/names.dmp 0 > ./Bacteria-Taxa.tsv
    rm taxdump.tar.gz

}


download_ftp $1
organize_files
